#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
    library(dplyr)
    library(pheatmap)
    library(survival)
    library(survminer)
    library(ggrepel)
    library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
cancer_type <- if (length(args) >= 1) args[1] else "TCGA-LUAD"

fig_dir <- file.path("results", "figures", "publication")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

PALETTE <- list(
    activated = "#E74C3C",
    repressed = "#2980B9",
    neutral   = "#95A5A6",
    highlight = "#F39C12",
    purple    = "#8E44AD",
    green     = "#2ECC71"
)

THEME_PUB <- theme_bw(base_size = 11) +
    theme(
        strip.background = element_rect(fill = "grey95"),
        panel.grid.minor = element_blank(),
        legend.key.size  = unit(0.4, "cm")
    )

save_pdf <- function(expr, path, width = 8, height = 6) {
    pdf(NULL)
    pdf(path, width = width, height = height)
    on.exit(dev.off())
    force(expr)
}

read_tsv_safe <- function(path) {
    if (!file.exists(path)) {
        return(NULL)
    }
    read.table(path,
        sep = "\t", header = TRUE, stringsAsFactors = FALSE,
        quote = "", fill = TRUE, comment.char = ""
    )
}

get_sym_col <- function(df) {
    found <- intersect(c("SYMBOL", "geneSymbol", "symbol"), colnames(df))
    if (length(found) == 0L) NULL else found[1]
}

integration_root <- file.path("pipeline_outputs", "integration")
run_dirs <- list.dirs(integration_root, recursive = FALSE, full.names = TRUE)

targets_long <- do.call(rbind, lapply(run_dirs, function(d) {
    f <- file.path(d, "direct_targets_all.tsv")
    if (!file.exists(f)) {
        return(NULL)
    }
    df <- read_tsv_safe(f)
    sym_col <- get_sym_col(df)
    gene_id_col <- if ("gene_id" %in% colnames(df)) "gene_id" else colnames(df)[1]
    data.frame(
        run = basename(d),
        symbol = if (!is.null(sym_col)) df[[sym_col]] else df[[gene_id_col]],
        log2FC = df$log2FoldChange,
        padj = df$padj,
        direction = ifelse(df$log2FoldChange > 0, "activated", "repressed"),
        stringsAsFactors = FALSE
    )
}))

conservation_df <- read_tsv_safe(
    file.path("pipeline_outputs", "multilinhagem", "target_conservation_matrix.tsv")
)
cox_df <- read_tsv_safe(file.path("results", "tables", "survival", cancer_type, "cox_results.tsv"))
immune_df <- read_tsv_safe(file.path("results", "tables", "immune_targets", "target_immune_classification.tsv"))
cat_enrich <- read_tsv_safe(file.path("results", "tables", "immune_targets", "category_enrichment.tsv"))

known_p53 <- c(
    "CDKN1A", "BAX", "MDM2", "BBC3", "GADD45A",
    "TIGAR", "GDF15", "PUMA", "TP53I3", "SESN1"
)

message("[fig1] Overview: targets per run and conservation")

if (!is.null(targets_long)) {
    p1a <- targets_long %>%
        count(run, direction) %>%
        ggplot(aes(x = run, y = n, fill = direction)) +
        geom_col(position = "stack", width = 0.6) +
        scale_fill_manual(values = c(
            activated = PALETTE$activated,
            repressed = PALETTE$repressed
        )) +
        labs(
            title = "A  Direct p53 targets per dataset",
            x = NULL, y = "Genes", fill = NULL
        ) +
        THEME_PUB +
        theme(
            axis.text.x = element_text(angle = 30, hjust = 1),
            legend.position = "top"
        )

    if (!is.null(conservation_df) && "n_runs" %in% colnames(conservation_df)) {
        p1b <- ggplot(conservation_df, aes(x = n_runs)) +
            geom_histogram(
                binwidth = 1, fill = PALETTE$green, color = "white",
                boundary = 0.5
            ) +
            scale_x_continuous(breaks = seq_len(max(conservation_df$n_runs, na.rm = TRUE))) +
            labs(
                title = "B  Target conservation across runs",
                x = "Number of runs", y = "Genes"
            ) +
            THEME_PUB

        save_pdf(
            print(p1a + p1b + plot_layout(widths = c(1.4, 1))),
            file.path(fig_dir, "fig1_overview.pdf"),
            width = 12, height = 5
        )
    } else {
        save_pdf(print(p1a), file.path(fig_dir, "fig1_overview.pdf"), width = 7, height = 5)
    }
}

message("[fig2] Volcano: all integration runs")

if (!is.null(targets_long)) {
    plot_df <- targets_long %>%
        filter(!is.na(padj) & !is.na(log2FC)) %>%
        mutate(
            neg_log10_padj = pmin(-log10(padj), 30),
            label = ifelse(symbol %in% known_p53, symbol, NA_character_),
            color_group = case_when(
                symbol %in% known_p53 ~ "Known p53 target",
                direction == "activated" ~ "Activated",
                direction == "repressed" ~ "Repressed",
                TRUE ~ "Other"
            )
        )

    color_vals <- c(
        "Known p53 target" = PALETTE$highlight,
        "Activated"        = PALETTE$activated,
        "Repressed"        = PALETTE$repressed,
        "Other"            = PALETTE$neutral
    )

    p2 <- ggplot(plot_df, aes(x = log2FC, y = neg_log10_padj, color = color_group)) +
        geom_point(alpha = 0.5, size = 1.2, shape = 16) +
        geom_text_repel(aes(label = label),
            size = 2.8, color = "black",
            max.overlaps = 20, segment.size = 0.3
        ) +
        scale_color_manual(values = color_vals) +
        geom_hline(
            yintercept = -log10(0.05), linetype = "dashed",
            linewidth = 0.4, color = "grey50"
        ) +
        geom_vline(
            xintercept = c(-1, 1), linetype = "dashed",
            linewidth = 0.4, color = "grey50"
        ) +
        facet_wrap(~run, scales = "free_y") +
        labs(
            title = "Volcano plots — direct p53 targets",
            x = "log2 Fold Change", y = "-log10(FDR)", color = NULL
        ) +
        THEME_PUB +
        theme(legend.position = "top")

    save_pdf(
        print(p2),
        file.path(fig_dir, "fig2_volcano.pdf"),
        width = max(7, 4 * length(unique(plot_df$run))), height = 5
    )
}

message("[fig3] Survival forest plot")

if (!is.null(cox_df)) {
    n_show <- min(25L, nrow(cox_df))
    forest_df <- cox_df[seq_len(n_show), ] %>%
        mutate(
            gene = factor(gene, levels = rev(gene)),
            sig  = padj < 0.05 & !is.na(padj)
        )

    p3 <- ggplot(
        forest_df,
        aes(x = hr, y = gene, xmin = hr_lo95, xmax = hr_hi95, color = sig)
    ) +
        geom_pointrange(size = 0.5, linewidth = 0.7) +
        geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
        scale_color_manual(
            values = c(`TRUE` = PALETTE$activated, `FALSE` = PALETTE$neutral),
            labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "n.s.")
        ) +
        labs(
            title = sprintf("Cox regression — %s", cancer_type),
            x = "Hazard Ratio (95 % CI)", y = NULL, color = NULL
        ) +
        THEME_PUB +
        theme(legend.position = "top")

    save_pdf(
        print(p3),
        file.path(fig_dir, "fig3_survival_forest.pdf"),
        width = 7, height = max(4, 0.28 * n_show + 2)
    )
}

message("[fig4] Immune category enrichment")

if (!is.null(cat_enrich)) {
    p4 <- cat_enrich %>%
        arrange(fold_enrichment) %>%
        mutate(
            category = factor(category, levels = category),
            sig = pvalue < 0.05
        ) %>%
        ggplot(aes(x = fold_enrichment, y = category, fill = sig)) +
        geom_col(width = 0.6) +
        geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
        scale_fill_manual(
            values = c(`TRUE` = PALETTE$activated, `FALSE` = PALETTE$neutral),
            labels = c(`TRUE` = "p < 0.05", `FALSE` = "n.s.")
        ) +
        labs(
            title = "Immune category enrichment of p53 targets",
            x = "Fold Enrichment", y = NULL, fill = NULL
        ) +
        THEME_PUB +
        theme(legend.position = "top")

    save_pdf(print(p4), file.path(fig_dir, "fig4_immune_enrichment.pdf"), width = 7, height = 5)
}

message("[fig5] Heatmap: conserved targets — direction across runs")

if (!is.null(conservation_df) && !is.null(targets_long)) {
    run_cols <- setdiff(
        colnames(conservation_df),
        c("symbol", "n_runs", "known_p53")
    )

    conserved_syms <- conservation_df$symbol[
        conservation_df$n_runs == length(run_cols) &
            !is.na(conservation_df$symbol)
    ]

    if (length(conserved_syms) >= 2L && length(run_cols) >= 2L) {
        lfc_mat <- matrix(NA_real_,
            nrow = length(conserved_syms), ncol = length(run_cols),
            dimnames = list(conserved_syms, run_cols)
        )

        for (run in run_cols) {
            sub <- targets_long[targets_long$run == run &
                targets_long$symbol %in% conserved_syms, ]
            lfc_mat[sub$symbol, run] <- sub$log2FC
        }

        keep <- names(sort(rowMeans(abs(lfc_mat), na.rm = TRUE),
            decreasing = TRUE
        ))[seq_len(min(40L, nrow(lfc_mat)))]
        lfc_mat <- lfc_mat[keep, , drop = FALSE]

        row_anno <- data.frame(
            Known_p53 = ifelse(rownames(lfc_mat) %in% known_p53, "yes", "no"),
            row.names = rownames(lfc_mat)
        )

        save_pdf(
            pheatmap(
                lfc_mat,
                annotation_row    = row_anno,
                annotation_colors = list(Known_p53 = c(yes = PALETTE$highlight, no = "grey85")),
                color             = colorRampPalette(c(PALETTE$repressed, "white", PALETTE$activated))(101),
                breaks            = seq(-3, 3, length.out = 102),
                border_color      = NA,
                fontsize_row      = 8,
                main              = "Conserved p53 targets — log2FC",
                cluster_cols      = TRUE,
                cluster_rows      = TRUE,
                silent            = TRUE
            ),
            file.path(fig_dir, "fig5_heatmap_conserved.pdf"),
            width = max(5, length(run_cols) * 2 + 2),
            height = max(6, nrow(lfc_mat) * 0.18 + 3)
        )
    }
}

message(sprintf("[done] Publication figures → %s", fig_dir))
