#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(survival)
    library(survminer)
    library(dplyr)
    library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

cancer_type <- if (length(args) >= 1) args[1] else "TCGA-BRCA"
max_genes_km <- if (length(args) >= 2) as.integer(args[2]) else 20L
tcga_dir <- if (length(args) >= 3) args[3] else file.path("data", "tcga")

out_dir <- file.path("results", "tables", "survival", cancer_type)
fig_dir <- file.path("results", "figures", "survival", cancer_type)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

message(sprintf("[setup] Cancer type : %s", cancer_type))
message(sprintf("[setup] TCGA dir    : %s", tcga_dir))

save_pdf <- function(expr, path, width = 8, height = 6) {
    pdf(path, width = width, height = height)
    on.exit(dev.off())
    force(expr)
}

read_tsv_safe <- function(path) {
    read.table(path,
        sep = "\t", header = TRUE, stringsAsFactors = FALSE,
        quote = "", fill = TRUE, comment.char = ""
    )
}

conserved_file <- file.path("pipeline_outputs", "multilinhagem", "conserved_targets_all_runs.tsv")
if (!file.exists(conserved_file)) {
    stop("File not found: ", conserved_file, "\nRun 10_multilinhagem.R first.", call. = FALSE)
}

conserved_df <- read_tsv_safe(conserved_file)
target_genes <- unique(conserved_df$symbol[!is.na(conserved_df$symbol) & nzchar(conserved_df$symbol)])
message(sprintf("[targets] Conserved p53 targets: %d", length(target_genes)))

manifest_path <- file.path(tcga_dir, "manifest.json")
clinical_path <- file.path(tcga_dir, "clinical.json")
expr_dir <- file.path(tcga_dir, "expression")

for (p in c(manifest_path, clinical_path, expr_dir)) {
    if (!file.exists(p)) {
        stop("Not found: ", p, "\nRun: python3 src/tcga_downloader.py --project ",
            cancer_type,
            call. = FALSE
        )
    }
}

message("[load] Reading manifest")
manifest <- jsonlite::fromJSON(manifest_path)
expr_files <- unlist(manifest)

missing <- expr_files[!file.exists(expr_files)]
if (length(missing) > 0L) {
    stop(length(missing), " expression file(s) listed in manifest not found on disk.\n",
        "Re-run: python3 src/tcga_downloader.py",
        call. = FALSE
    )
}

message(sprintf("[load] Reading %d expression file(s)", length(expr_files)))

read_expr_file <- function(path) {
    df <- read.table(path,
        sep = "\t", header = TRUE, stringsAsFactors = FALSE,
        comment.char = "#", fill = TRUE, quote = ""
    )
    df <- df[!grepl("^N_", df[[1]]), ]
    df
}

first_df <- read_expr_file(expr_files[1])
gene_col <- "gene_name"
count_col <- "tpm_unstranded"

if (!gene_col %in% colnames(first_df)) stop("Column 'gene_name' not found in expression files.", call. = FALSE)
if (!count_col %in% colnames(first_df)) stop("Column 'tpm_unstranded' not found in expression files.", call. = FALSE)

gene_names <- first_df[[gene_col]]
target_idx <- which(gene_names %in% target_genes)

if (length(target_idx) == 0L) {
    stop("None of the conserved target genes found in expression files.", call. = FALSE)
}

found_genes <- gene_names[target_idx]
message(sprintf(
    "[match] Targets found in expression matrix: %d / %d",
    length(unique(found_genes)), length(target_genes)
))

sample_ids <- basename(expr_files)
sample_ids <- sub("\\.rna_seq.*$", "", sample_ids)

expr_mat <- matrix(NA_real_,
    nrow = length(target_idx), ncol = length(expr_files),
    dimnames = list(found_genes, sample_ids)
)

for (i in seq_along(expr_files)) {
    df <- read_expr_file(expr_files[i])
    df <- df[!grepl("^N_", df[[gene_col]]), ]
    idx <- match(target_genes, df[[gene_col]])
    vals <- as.numeric(df[[count_col]][idx[!is.na(idx)]])
    expr_mat[!is.na(idx), i] <- vals
}

message("[load] Reading clinical data")
clin_list <- jsonlite::fromJSON(clinical_path, simplifyDataFrame = FALSE)

clin_df <- do.call(rbind, lapply(clin_list, function(x) {
    data.frame(
        submitter_id           = x$submitter_id %||% NA_character_,
        vital_status           = x$vital_status %||% NA_character_,
        days_to_death          = x$days_to_death %||% NA_real_,
        days_to_last_follow_up = x$days_to_last_follow_up %||% NA_real_,
        age_at_diagnosis       = x$age_at_diagnosis %||% NA_real_,
        stringsAsFactors       = FALSE
    )
}))

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0L && !is.na(a)) a else b

clin_df <- do.call(rbind, lapply(clin_list, function(x) {
    data.frame(
        submitter_id           = if (!is.null(x$submitter_id)) x$submitter_id else NA_character_,
        vital_status           = if (!is.null(x$vital_status)) x$vital_status else NA_character_,
        days_to_death          = if (!is.null(x$days_to_death)) as.numeric(x$days_to_death) else NA_real_,
        days_to_last_follow_up = if (!is.null(x$days_to_last_follow_up)) as.numeric(x$days_to_last_follow_up) else NA_real_,
        age_at_diagnosis       = if (!is.null(x$age_at_diagnosis)) as.numeric(x$age_at_diagnosis) else NA_real_,
        stringsAsFactors       = FALSE
    )
}))

message(sprintf("[clinical] %d records loaded", nrow(clin_df)))

clin_df$os_time <- ifelse(
    !is.na(clin_df$days_to_death),
    clin_df$days_to_death,
    clin_df$days_to_last_follow_up
)
clin_df$os_event <- as.integer(toupper(clin_df$vital_status) == "DEAD")

common_ids <- intersect(sample_ids, clin_df$submitter_id)
message(sprintf("[merge] Samples with expression + clinical data: %d", length(common_ids)))

if (length(common_ids) < 10L) {
    stop("Too few matched samples (", length(common_ids),
        "). The manifest file IDs may not match submitter_ids in clinical.json.\n",
        "Re-run: python3 src/tcga_downloader.py --project ", cancer_type,
        call. = FALSE
    )
}

expr_col_idx <- match(common_ids, sample_ids)
clin_row_idx <- match(common_ids, clin_df$submitter_id)

analysis_df <- as.data.frame(t(expr_mat[, expr_col_idx, drop = FALSE]))
rownames(analysis_df) <- common_ids
analysis_df$os_time <- clin_df$os_time[clin_row_idx]
analysis_df$os_event <- clin_df$os_event[clin_row_idx]

analysis_df <- analysis_df[!is.na(analysis_df$os_time) & analysis_df$os_time > 0, ]
message(sprintf("[survival] Samples with valid OS data: %d", nrow(analysis_df)))

message("[cox] Running univariate Cox regression")

unique_genes <- unique(found_genes)

cox_results <- do.call(rbind, lapply(unique_genes, function(gene) {
    if (!gene %in% colnames(analysis_df)) {
        return(NULL)
    }

    expr <- analysis_df[[gene]]
    if (all(is.na(expr)) || var(expr, na.rm = TRUE) == 0) {
        return(NULL)
    }

    fit <- tryCatch(
        coxph(Surv(os_time, os_event) ~ expr, data = analysis_df),
        error = function(e) NULL
    )
    if (is.null(fit)) {
        return(NULL)
    }

    s <- summary(fit)
    data.frame(
        gene = gene,
        hr = s$coefficients[1, "exp(coef)"],
        hr_lo95 = s$conf.int[1, "lower .95"],
        hr_hi95 = s$conf.int[1, "upper .95"],
        z = s$coefficients[1, "z"],
        pvalue = s$coefficients[1, "Pr(>|z|)"],
        stringsAsFactors = FALSE
    )
}))

if (is.null(cox_results) || nrow(cox_results) == 0L) {
    stop("Cox regression produced no results. Check expression and clinical data.", call. = FALSE)
}

cox_results$padj <- p.adjust(cox_results$pvalue, method = "BH")
cox_results <- cox_results[order(cox_results$pvalue), ]

message(sprintf("[cox] Genes with FDR < 0.05: %d", sum(cox_results$padj < 0.05, na.rm = TRUE)))

write.table(cox_results, file.path(out_dir, "cox_results.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)

message("[plots] Forest plot")

n_forest <- min(30L, nrow(cox_results))
forest_df <- cox_results[seq_len(n_forest), ]
forest_df$gene <- factor(forest_df$gene, levels = rev(forest_df$gene))

save_pdf(
    print(
        ggplot(forest_df, aes(
            x = hr, y = gene, xmin = hr_lo95, xmax = hr_hi95,
            color = padj < 0.05
        )) +
            geom_pointrange(size = 0.5) +
            geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
            scale_color_manual(
                values = c(`TRUE` = "#E74C3C", `FALSE` = "grey40"),
                labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "FDR >= 0.05")
            ) +
            labs(
                title = sprintf("Univariate Cox — %s (top %d targets)", cancer_type, n_forest),
                x = "Hazard Ratio (95% CI)", y = NULL, color = NULL
            ) +
            theme_bw(base_size = 11) +
            theme(legend.position = "top")
    ),
    file.path(fig_dir, "forest_plot.pdf"),
    width = 7, height = max(4, 0.25 * n_forest + 2)
)

sig_genes <- cox_results$gene[cox_results$padj < 0.05 & !is.na(cox_results$padj)]
km_genes <- sig_genes[seq_len(min(max_genes_km, length(sig_genes)))]

if (length(km_genes) == 0L) {
    message("[plots] No FDR-significant genes — skipping KM curves")
} else {
    message(sprintf("[plots] KM curves for %d significant gene(s)", length(km_genes)))

    for (gene in km_genes) {
        expr <- analysis_df[[gene]]
        km_df <- data.frame(
            os_time  = analysis_df$os_time,
            os_event = analysis_df$os_event,
            group    = ifelse(expr >= median(expr, na.rm = TRUE), "High", "Low")
        )

        fit <- survfit(Surv(os_time, os_event) ~ group, data = km_df)
        hr_row <- cox_results[cox_results$gene == gene, ]
        subtitle <- sprintf(
            "HR = %.2f (%.2f-%.2f), p = %.3g, FDR = %.3g",
            hr_row$hr, hr_row$hr_lo95, hr_row$hr_hi95,
            hr_row$pvalue, hr_row$padj
        )

        p <- ggsurvplot(
            fit,
            data = km_df, pval = FALSE, risk.table = TRUE,
            palette = c(High = "#E74C3C", Low = "#2980B9"),
            title = sprintf("%s expression — %s", gene, cancer_type),
            subtitle = subtitle,
            xlab = "Days", ylab = "Overall survival probability",
            legend.title = gene,
            ggtheme = theme_bw(base_size = 11)
        )

        safe_gene <- gsub("[^A-Za-z0-9_]", "_", gene)
        save_pdf(print(p), file.path(fig_dir, sprintf("km_%s.pdf", safe_gene)),
            width = 7, height = 6
        )
    }
}

message(sprintf("[done] Cox results -> %s", out_dir))
message(sprintf("[done] Figures     -> %s", fig_dir))
