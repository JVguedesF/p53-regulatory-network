#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

conserved_path  <- if (length(args) >= 1) args[1] else
    file.path("pipeline_outputs", "multilinhagem", "conserved_targets_all_runs.tsv")
survival_path   <- if (length(args) >= 2 && nzchar(args[2])) args[2] else NULL
cancer_type     <- if (length(args) >= 3) args[3] else "TCGA-BRCA"
fdr_survival    <- if (length(args) >= 4) as.numeric(args[4]) else 0.05
background_size <- if (length(args) >= 5) as.integer(args[5]) else 20000L

out_dir <- file.path("results", "tables",  "immune_targets")
fig_dir <- file.path("results", "figures", "immune_targets")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

save_pdf <- function(expr, path, width = 8, height = 6) {
    pdf(path, width = width, height = height)
    on.exit(dev.off())
    force(expr)
}

IMMUNE_SETS <- list(
    checkpoint          = c("PDCD1","CD274","CTLA4","LAG3","HAVCR2","TIGIT","VSIR",
                             "PDCD1LG2","CD80","CD86","ICOS","ICOSLG"),
    cytokines_receptors = c("IFNG","TNF","IL2","IL6","IL10","IL12A","IL12B",
                             "IL15","IL18","TGFB1","CXCL9","CXCL10","CXCL11",
                             "CCL2","CCL5","IFNGR1","IFNGR2","IL2RA"),
    antigen_presentation = c("HLA-A","HLA-B","HLA-C","HLA-DRA","HLA-DRB1",
                              "B2M","TAP1","TAP2","TAPBP","NLRC5"),
    innate_sensing      = c("STING1","CGAS","IRF3","IRF7","MAVS","DDX58",
                             "IFIH1","TLR3","TLR4","TLR9","MYD88"),
    nk_cytotoxicity     = c("KLRK1","KLRD1","NCR1","NCR3","MICA","MICB",
                             "ULBP1","ULBP2","ULBP3","RAET1E"),
    t_cell_trafficking  = c("CXCR3","CCR5","CXCR6","CCR2","SELL",
                             "ITGAL","ITGB2","ICAM1","VCAM1"),
    apoptosis_immune    = c("FASLG","FAS","TRAIL","TNFRSF10A","TNFRSF10B",
                             "CASP8","CASP3","CFLAR")
)

all_immune_genes <- unique(unlist(IMMUNE_SETS))

if (!file.exists(conserved_path)) {
    stop("File not found: ", conserved_path, "\nRun 10_multilinhagem.R first.", call. = FALSE)
}

conserved_df <- read.table(conserved_path, sep = "\t", header = TRUE,
                            stringsAsFactors = FALSE, quote = "", fill = TRUE)
target_genes <- unique(conserved_df$symbol[!is.na(conserved_df$symbol) & nzchar(conserved_df$symbol)])
message(sprintf("[targets] Conserved p53 targets: %d", length(target_genes)))
message(sprintf("[immune]  Immune gene database : %d", length(all_immune_genes)))

survival_genes <- NULL
if (!is.null(survival_path) && file.exists(survival_path)) {
    cox_df <- read.table(survival_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    survival_genes <- cox_df$gene[!is.na(cox_df$padj) & cox_df$padj < fdr_survival]
    message(sprintf("[survival] Survival-significant genes (FDR < %.2f): %d",
                    fdr_survival, length(survival_genes)))
}

classified_df <- do.call(rbind, lapply(target_genes, function(gene) {
    matched <- names(IMMUNE_SETS)[sapply(IMMUNE_SETS, function(m) gene %in% m)]
    data.frame(
        gene              = gene,
        is_immune         = length(matched) > 0L,
        immune_categories = paste(matched, collapse = ";"),
        stringsAsFactors  = FALSE
    )
}))

immune_targets <- classified_df$gene[classified_df$is_immune]
message(sprintf("[overlap] Immune targets: %d / %d", length(immune_targets), length(target_genes)))

write.table(classified_df, file.path(out_dir, "target_immune_classification.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

fisher_enrichment <- function(targets, immune_bg, bg_size) {
    k <- length(intersect(targets, immune_bg))
    n <- length(targets)
    K <- length(immune_bg)
    N <- bg_size
    mat <- matrix(c(k, n - k, K - k, max(N - K - (n - k), 0L)), nrow = 2)
    p   <- fisher.test(mat, alternative = "greater")$p.value
    list(overlap = k, n_targets = n, n_immune = K, background = N,
         expected = round(n * K / N, 2),
         fold_enrichment = if (K > 0 && n > 0) round((k / n) / (K / N), 3) else 0,
         pvalue_fisher = round(p, 6))
}

global_enr <- fisher_enrichment(target_genes, all_immune_genes, background_size)
message(sprintf("[enrichment] Overlap=%d, FE=%.3f, p=%.4g",
                global_enr$overlap, global_enr$fold_enrichment, global_enr$pvalue_fisher))

cat_df <- do.call(rbind, lapply(names(IMMUNE_SETS), function(cat) {
    members <- IMMUNE_SETS[[cat]]
    res     <- fisher_enrichment(target_genes, members, background_size)
    data.frame(category = cat, n_category = length(members),
               overlap = res$overlap, genes = paste(intersect(target_genes, members), collapse = ";"),
               fold_enrichment = res$fold_enrichment, pvalue = res$pvalue_fisher,
               stringsAsFactors = FALSE)
}))
cat_df <- cat_df[order(cat_df$pvalue), ]
cat_df$padj <- p.adjust(cat_df$pvalue, method = "BH")

write.table(cat_df, file.path(out_dir, "category_enrichment.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

if (!is.null(survival_genes)) {
    both        <- intersect(immune_targets, survival_genes)
    priority_df <- classified_df[classified_df$gene %in% both, ]
    priority_df$survival_significant <- TRUE
    write.table(priority_df, file.path(out_dir, "priority_immune_survival_targets.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    message(sprintf("[priority] Immune + survival-significant: %d", length(both)))
}

message("[plots] Category enrichment bar chart")

plot_df <- cat_df %>%
    mutate(category = factor(category, levels = rev(category)))

save_pdf(
    print(
        ggplot(plot_df, aes(x = fold_enrichment, y = category, fill = pvalue < 0.05)) +
            geom_col(width = 0.6) +
            geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
            scale_fill_manual(values = c(`TRUE` = "#E74C3C", `FALSE` = "#95A5A6"),
                              labels = c(`TRUE` = "p < 0.05", `FALSE` = "n.s.")) +
            labs(title = "Immune Category Enrichment of p53 Targets",
                 x = "Fold Enrichment", y = NULL, fill = NULL) +
            theme_bw(base_size = 11) +
            theme(legend.position = "top")
    ),
    file.path(fig_dir, "category_enrichment.pdf"),
    width = 7, height = max(3, 0.5 * nrow(cat_df) + 1)
)

message("[plots] Immune target heatmap")

immune_overlap <- intersect(target_genes, all_immune_genes)

if (length(immune_overlap) > 0L) {
    mat_df <- do.call(rbind, lapply(immune_overlap, function(gene) {
        sapply(names(IMMUNE_SETS), function(cat) as.integer(gene %in% IMMUNE_SETS[[cat]]))
    }))
    rownames(mat_df) <- immune_overlap

    heat_df <- as.data.frame(mat_df) %>%
        tibble::rownames_to_column("gene") %>%
        tidyr::pivot_longer(-gene, names_to = "category", values_to = "member")

    save_pdf(
        print(
            ggplot(heat_df, aes(x = category, y = gene, fill = factor(member))) +
                geom_tile(color = "white", linewidth = 0.3) +
                scale_fill_manual(values = c("0" = "grey92", "1" = "#2980B9"), guide = "none") +
                labs(title = "p53 Immune Targets x Category Membership",
                     x = NULL, y = NULL) +
                theme_bw(base_size = 10) +
                theme(axis.text.x = element_text(angle = 35, hjust = 1))
        ),
        file.path(fig_dir, "immune_target_heatmap.pdf"),
        width  = max(5, length(IMMUNE_SETS) * 0.9),
        height = max(4, length(immune_overlap) * 0.35 + 1)
    )
}

if (!is.null(survival_genes)) {
    message("[plots] Immune vs survival overlap")

    counts_df <- data.frame(
        group = c("Immune Targets", "Survival Sig."),
        n     = c(length(immune_targets), length(survival_genes))
    )
    both_genes <- intersect(immune_targets, survival_genes)

    save_pdf(
        {
            p <- ggplot(counts_df, aes(x = group, y = n, fill = group)) +
                geom_col(width = 0.5, show.legend = FALSE) +
                scale_fill_manual(values = c("Immune Targets" = "#2ECC71",
                                             "Survival Sig."  = "#E74C3C")) +
                labs(title = "Immune Targets with Survival Significance",
                     x = NULL, y = "Gene count") +
                theme_bw(base_size = 11)
            if (length(both_genes) > 0L) {
                p <- p + annotate("text", x = 1.5, y = max(counts_df$n) * 0.5,
                                  label = paste("Overlap:\n", paste(sort(both_genes), collapse = "\n")),
                                  size = 3.5, hjust = 0.5,
                                  bbox = list(boxstyle = "round", facecolor = "#FAD7A0"))
            }
            print(p)
        },
        file.path(fig_dir, "immune_survival_overlap.pdf"),
        width = 7, height = 5
    )
}

message(sprintf("[done] Tables  -> %s", out_dir))
message(sprintf("[done] Figures -> %s", fig_dir))