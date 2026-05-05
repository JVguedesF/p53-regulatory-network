#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(TCGAbiolinks)
    library(SummarizedExperiment)
    library(survival)
    library(survminer)
    library(dplyr)
    library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

cancer_type <- if (length(args) >= 1) args[1] else "TCGA-LUAD"
max_genes_km <- if (length(args) >= 2) as.integer(args[2]) else 20L
cache_dir <- if (length(args) >= 3) args[3] else file.path("data", "tcga_cache")

out_dir <- file.path("results", "survival", cancer_type)
fig_dir <- file.path("results", "figures", "survival", cancer_type)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

message(sprintf("[setup] Cancer type : %s", cancer_type))
message(sprintf("[setup] Cache dir   : %s", cache_dir))

conserved_file <- file.path("results", "multilinhagem", "conserved_targets_all_runs.tsv")

if (!file.exists(conserved_file)) {
    stop(
        "File not found: ", conserved_file,
        "\nRun 10_multilinhagem.R first.",
        call. = FALSE
    )
}

conserved_df <- read.table(conserved_file,
    sep = "\t", header = TRUE, stringsAsFactors = FALSE,
    quote = "", fill = TRUE, comment.char = ""
)

if (!"symbol" %in% colnames(conserved_df)) {
    stop("Expected column 'symbol' not found in ", conserved_file, call. = FALSE)
}

target_genes <- unique(conserved_df$symbol[!is.na(conserved_df$symbol)])
message(sprintf("[targets] Conserved p53 targets to test: %d", length(target_genes)))

rds_path <- file.path(cache_dir, sprintf("%s_rnaseq.rds", cancer_type))

if (file.exists(rds_path)) {
    message("[tcga] Loading cached expression data")
    se <- readRDS(rds_path)
} else {
    message("[tcga] Querying GDC for RNA-seq data (this may take a few minutes)")

    query <- GDCquery(
        project           = cancer_type,
        data.category     = "Transcriptome Profiling",
        data.type         = "Gene Expression Quantification",
        workflow.type     = "STAR - Counts",
        sample.type       = "Primary Tumor"
    )

    GDCdownload(query, directory = cache_dir)
    se <- GDCprepare(query, directory = cache_dir)
    saveRDS(se, rds_path)
    message(sprintf("[tcga] Data cached to %s", rds_path))
}

message(sprintf("[tcga] %d samples, %d genes", ncol(se), nrow(se)))

clin_rds <- file.path(cache_dir, sprintf("%s_clinical.rds", cancer_type))

if (file.exists(clin_rds)) {
    message("[tcga] Loading cached clinical data")
    clin_df <- readRDS(clin_rds)
} else {
    message("[tcga] Downloading clinical data")
    clin_df <- GDCquery_clinic(project = cancer_type, type = "clinical")
    saveRDS(clin_df, clin_rds)
}

os_time_candidates <- c("days_to_death", "days_to_last_follow_up", "overall_survival")
os_status_candidates <- c("vital_status", "overall_survival_status")

os_time_col <- intersect(os_time_candidates, colnames(clin_df))[1]
os_status_col <- intersect(os_status_candidates, colnames(clin_df))[1]

if (is.na(os_time_col) || is.na(os_status_col)) {
    stop(
        "Could not identify OS time/status columns in clinical data.\n",
        "Available columns: ", paste(colnames(clin_df), collapse = ", "),
        call. = FALSE
    )
}

message(sprintf("[clinical] Using '%s' (time) and '%s' (status)", os_time_col, os_status_col))

se_barcodes <- substr(colnames(se), 1, 12)
clin_barcodes <- clin_df$submitter_id

common <- intersect(se_barcodes, clin_barcodes)
message(sprintf("[merge] Samples with both expression and clinical data: %d", length(common)))

if (length(common) < 10L) {
    stop("Too few matched samples (", length(common), "). Check TCGA project code.", call. = FALSE)
}

se_idx <- which(se_barcodes %in% common)
clin_idx <- match(se_barcodes[se_idx], clin_barcodes)

assay_name <- if ("tpm_unstrand" %in% assayNames(se)) "tpm_unstrand" else assayNames(se)[1]
message(sprintf("[expression] Using assay: %s", assay_name))

expr_mat <- assay(se, assay_name)[, se_idx]

gene_meta <- rowData(se)
sym_col_se <- if ("gene_name" %in% colnames(gene_meta)) "gene_name" else NULL

if (is.null(sym_col_se)) {
    stop("Could not find 'gene_name' in rowData. Inspect rowData(se) columns.", call. = FALSE)
}

gene_symbols_se <- gene_meta[[sym_col_se]]
target_idx <- which(gene_symbols_se %in% target_genes)

if (length(target_idx) == 0L) {
    stop("None of the conserved target genes were found in the TCGA expression matrix.", call. = FALSE)
}

found_genes <- gene_symbols_se[target_idx]
message(sprintf(
    "[match] Conserved targets found in TCGA matrix: %d / %d",
    length(unique(found_genes)), length(target_genes)
))

expr_targets <- t(expr_mat[target_idx, ])
colnames(expr_targets) <- found_genes

os_time <- as.numeric(clin_df[[os_time_col]][clin_idx])
os_status <- toupper(clin_df[[os_status_col]][clin_idx])
os_event <- as.integer(os_status %in% c("DEAD", "1", "TRUE"))

analysis_df <- data.frame(
    sample = colnames(se)[se_idx],
    os_time = os_time,
    os_event = os_event,
    expr_targets,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

analysis_df <- analysis_df[!is.na(analysis_df$os_time) & analysis_df$os_time > 0, ]
message(sprintf("[survival] Samples with valid OS data: %d", nrow(analysis_df)))

message("[cox] Running univariate Cox regression for each target gene")

unique_genes <- unique(found_genes)

cox_results <- do.call(rbind, lapply(unique_genes, function(gene) {
    cols <- which(colnames(analysis_df) == gene)
    expr <- if (length(cols) == 1L) {
        analysis_df[[cols]]
    } else {
        rowMeans(analysis_df[, cols], na.rm = TRUE)
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

cox_results$padj <- p.adjust(cox_results$pvalue, method = "BH")
cox_results <- cox_results[order(cox_results$pvalue), ]

message(sprintf("[cox] Genes with FDR < 0.05: %d", sum(cox_results$padj < 0.05, na.rm = TRUE)))

write.table(cox_results,
    file.path(out_dir, "cox_results.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)

message("[plots] Forest plot")

n_forest <- min(30L, nrow(cox_results))
forest_df <- cox_results[seq_len(n_forest), ]
forest_df$gene <- factor(forest_df$gene, levels = rev(forest_df$gene))

save_pdf <- function(expr, path, width = 8, height = 6) {
    pdf(NULL)
    pdf(path, width = width, height = height)
    on.exit(dev.off())
    force(expr)
}

save_pdf(
    print(
        ggplot(forest_df, aes(
            x = hr, y = gene,
            xmin = hr_lo95, xmax = hr_hi95,
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
                x     = "Hazard Ratio (95% CI)",
                y     = NULL,
                color = NULL
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
        cols <- which(colnames(analysis_df) == gene)
        expr <- if (length(cols) == 1L) {
            analysis_df[[cols]]
        } else {
            rowMeans(analysis_df[, cols], na.rm = TRUE)
        }

        km_df <- data.frame(
            os_time  = analysis_df$os_time,
            os_event = analysis_df$os_event,
            group    = ifelse(expr >= median(expr, na.rm = TRUE), "High", "Low")
        )

        fit <- survfit(Surv(os_time, os_event) ~ group, data = km_df)

        hr_row <- cox_results[cox_results$gene == gene, ]
        subtitle <- sprintf(
            "HR = %.2f (%.2f–%.2f), p = %.3g, FDR = %.3g",
            hr_row$hr, hr_row$hr_lo95, hr_row$hr_hi95,
            hr_row$pvalue, hr_row$padj
        )

        p <- ggsurvplot(
            fit,
            data          = km_df,
            pval          = FALSE,
            risk.table    = TRUE,
            palette       = c(High = "#E74C3C", Low = "#2980B9"),
            title         = sprintf("%s expression — %s", gene, cancer_type),
            subtitle      = subtitle,
            xlab          = "Days",
            ylab          = "Overall survival probability",
            legend.title  = sprintf("%s", gene),
            ggtheme       = theme_bw(base_size = 11)
        )

        safe_gene <- gsub("[^A-Za-z0-9_]", "_", gene)
        save_pdf(
            print(p),
            file.path(fig_dir, sprintf("km_%s.pdf", safe_gene)),
            width = 7, height = 6
        )
    }
}

message(sprintf("[done] Cox results → %s", out_dir))
message(sprintf("[done] Figures     → %s", fig_dir))
