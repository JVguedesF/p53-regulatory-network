#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(dplyr)
    library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    stop(
        "Usage: Rscript 10_integration.R",
        " <chipseq_tsv> <rnaseq_tsv>",
        " [chip_dataset_name] [rna_dataset_name]",
        " [lfc_threshold=1] [padj_threshold=0.05] [promoter_window=3000]",
        call. = FALSE
    )
}

chipseq_tsv <- args[1]
rnaseq_tsv <- args[2]
chip_base <- if (length(args) >= 3) args[3] else sub("\\.tsv$", "", basename(chipseq_tsv))
rna_base <- if (length(args) >= 4) args[4] else sub("\\.tsv$", "", basename(rnaseq_tsv))
lfc_threshold <- if (length(args) >= 5) as.numeric(args[5]) else 1.0
padj_threshold <- if (length(args) >= 6) as.numeric(args[6]) else 0.05
promoter_window <- if (length(args) >= 7) as.integer(args[7]) else 3000L

run_label <- sprintf("%s_%s", chip_base, rna_base)

# integration TSVs feed scripts 11, 12, 13 → pipeline_outputs
out_dir <- file.path("pipeline_outputs", "integration", run_label)
# network is a final deliverable (Cytoscape) → results
net_dir <- file.path("results", "networks")
# figures are for humans → results
fig_dir <- file.path("results", "figures", "integration", run_label)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(net_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

save_pdf <- function(expr, path, width = 8, height = 6) {
    pdf(NULL)
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

message(sprintf("[setup] ChIP-seq dataset : %s", chip_base))
message(sprintf("[setup] RNA-seq dataset  : %s", rna_base))

anno_dir <- file.path("pipeline_outputs", "chipseq", chip_base, "annotation")
if (!dir.exists(anno_dir)) {
    stop("Annotation directory not found: ", anno_dir,
        "\nRun script 07 first.",
        call. = FALSE
    )
}

anno_files <- list.files(anno_dir,
    pattern = "_annotation\\.tsv$",
    recursive = TRUE, full.names = TRUE
)
if (length(anno_files) == 0L) {
    stop("No annotation TSV files found in: ", anno_dir,
        "\nRun script 07 first.",
        call. = FALSE
    )
}

message(sprintf("[peaks] Loading %d annotation file(s)", length(anno_files)))

anno_df <- do.call(rbind, lapply(anno_files, function(f) {
    df <- read_tsv_safe(f)
    df$source_file <- basename(dirname(f))
    df
}))

message(sprintf("[peaks] Total annotated peaks: %d", nrow(anno_df)))

degs_dir <- file.path("pipeline_outputs", "rnaseq", rna_base, "deseq2")
if (!dir.exists(degs_dir)) {
    stop("DESeq2 results not found: ", degs_dir,
        "\nRun script 09 first.",
        call. = FALSE
    )
}

degs_files <- list.files(degs_dir,
    pattern = "degs_all\\.tsv$",
    recursive = TRUE, full.names = TRUE
)
if (length(degs_files) == 0L) {
    stop("No degs_all.tsv found in: ", degs_dir,
        "\nRun script 09 first.",
        call. = FALSE
    )
}

message(sprintf("[degs] Loading DESeq2 results from: %s", degs_files[1]))
degs_df <- read_tsv_safe(degs_files[1])
message(sprintf("[degs] Total tested genes: %d", nrow(degs_df)))

promoter_peaks <- anno_df[grepl("Promoter", anno_df$annotation, ignore.case = TRUE), ]
message(sprintf("[peaks] Peaks in promoters (±%dkb TSS): %d", promoter_window / 1000, nrow(promoter_peaks)))

degs_sig <- degs_df[
    !is.na(degs_df$padj) &
        degs_df$padj < padj_threshold &
        !is.na(degs_df$log2FoldChange) &
        abs(degs_df$log2FoldChange) > lfc_threshold,
]
message(sprintf(
    "[degs] Significant DEGs (|lfc|>%.1f, padj<%.2f): %d",
    lfc_threshold, padj_threshold, nrow(degs_sig)
))

peak_gene_col <- if ("geneId" %in% colnames(promoter_peaks)) {
    "geneId"
} else if ("ENTREZID" %in% colnames(promoter_peaks)) {
    "ENTREZID"
} else {
    stop("Cannot find gene ID column in annotation TSV.", call. = FALSE)
}

deg_gene_col <- if ("gene_id" %in% colnames(degs_df)) {
    "gene_id"
} else if ("ENSEMBL" %in% colnames(degs_df)) {
    "ENSEMBL"
} else {
    stop("Cannot find gene ID column in DESeq2 results.", call. = FALSE)
}

direct_targets <- merge(
    promoter_peaks, degs_sig,
    by.x = peak_gene_col, by.y = deg_gene_col,
    all = FALSE
)

indirect_targets <- degs_sig[
    !degs_sig[[deg_gene_col]] %in% promoter_peaks[[peak_gene_col]],
]

bound_no_change <- promoter_peaks[
    !promoter_peaks[[peak_gene_col]] %in% degs_sig[[deg_gene_col]],
]

direct_activated <- direct_targets[!is.na(direct_targets$log2FoldChange) & direct_targets$log2FoldChange > 0, ]
direct_repressed <- direct_targets[!is.na(direct_targets$log2FoldChange) & direct_targets$log2FoldChange < 0, ]

message(sprintf("[integration] Direct targets (peak + DEG): %d", nrow(direct_targets)))
message(sprintf("[integration]   Activated : %d", nrow(direct_activated)))
message(sprintf("[integration]   Repressed : %d", nrow(direct_repressed)))
message(sprintf("[integration] Indirect targets (DEG only): %d", nrow(indirect_targets)))
message(sprintf("[integration] Bound, no expression change : %d", nrow(bound_no_change)))

write_tsv <- function(df, path) {
    write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

write_tsv(direct_targets, file.path(out_dir, "direct_targets_all.tsv"))
write_tsv(direct_activated, file.path(out_dir, "direct_targets_activated.tsv"))
write_tsv(direct_repressed, file.path(out_dir, "direct_targets_repressed.tsv"))
write_tsv(indirect_targets, file.path(out_dir, "indirect_targets.tsv"))
write_tsv(bound_no_change, file.path(out_dir, "bound_no_change.tsv"))

write_tsv(
    data.frame(
        Category = c("Direct_Activated", "Direct_Repressed", "Indirect_DEG", "Bound_No_Change"),
        N = c(
            nrow(direct_activated), nrow(direct_repressed),
            nrow(indirect_targets), nrow(bound_no_change)
        )
    ),
    file.path(out_dir, "integration_summary.tsv")
)

score_col <- if ("V5" %in% colnames(direct_targets)) "V5" else if ("score" %in% colnames(direct_targets)) "score" else NULL

message("[plots] Scatter: peak score vs log2FC")

if (!is.null(score_col) && nrow(direct_targets) > 1L) {
    cor_res <- tryCatch(
        cor.test(direct_targets[[score_col]], abs(direct_targets$log2FoldChange), method = "spearman"),
        error = function(e) NULL
    )
    cor_label <- if (!is.null(cor_res)) {
        sprintf("Spearman r = %.3f, p = %.2e", cor_res$estimate, cor_res$p.value)
    } else {
        ""
    }

    plot_df <- direct_targets
    plot_df$direction <- ifelse(plot_df$log2FoldChange > 0, "Activated", "Repressed")
    gene_sym_col <- if ("SYMBOL" %in% colnames(plot_df)) "SYMBOL" else if ("geneSymbol" %in% colnames(plot_df)) "geneSymbol" else NULL
    key_genes <- c("CDKN1A", "BAX", "MDM2", "BBC3", "GADD45A", "TIGAR")

    save_pdf(
        {
            p <- ggplot(plot_df, aes(x = .data[[score_col]], y = log2FoldChange, color = direction)) +
                geom_point(alpha = 0.6, size = 1.5, shape = 16) +
                scale_color_manual(values = c(Activated = "#E74C3C", Repressed = "#2980B9")) +
                geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey50") +
                labs(
                    title = sprintf("Peak Intensity vs Transcriptional Response\n%s", cor_label),
                    x     = "Peak Score (MACS2)",
                    y     = "Log2 Fold Change (RNA-seq)",
                    color = NULL
                ) +
                theme_bw(base_size = 12) +
                theme(legend.position = "top")

            if (!is.null(gene_sym_col)) {
                label_df <- plot_df[!is.na(plot_df[[gene_sym_col]]) & plot_df[[gene_sym_col]] %in% key_genes, ]
                if (nrow(label_df) > 0L) {
                    p <- p + ggrepel::geom_text_repel(
                        data = label_df, aes(label = .data[[gene_sym_col]]),
                        size = 3, color = "black", max.overlaps = 20
                    )
                }
            }
            print(p)
        },
        file.path(fig_dir, "peak_score_vs_log2fc.pdf")
    )
}

message("[plots] Venn: promoter peaks ∩ DEGs")

tryCatch(
    {
        library(ggVennDiagram, warn.conflicts = FALSE)
        venn_list <- list(
            Promoter_Peaks = as.character(promoter_peaks[[peak_gene_col]]),
            DEGs           = as.character(degs_sig[[deg_gene_col]])
        )
        save_pdf(
            print(ggVennDiagram(venn_list,
                label_alpha = 0,
                category.names = c("Promoter Peaks", "DEGs")
            ) +
                scale_fill_gradient(low = "#FFFFFF", high = "#2980B9") +
                ggtitle("ChIP-seq Promoter Peaks ∩ RNA-seq DEGs")),
            file.path(fig_dir, "venn_peaks_degs.pdf"),
            width = 6, height = 6
        )
    },
    error = function(e) message("[warn] Venn diagram skipped: ", conditionMessage(e))
)

message("[network] Building p53 regulatory network")

tryCatch(
    {
        library(igraph, warn.conflicts = FALSE)

        edges_act <- data.frame(
            from = chip_base,
            to = as.character(direct_activated[[peak_gene_col]]),
            type = "activation",
            weight = if (!is.null(score_col)) direct_activated[[score_col]] else 1,
            log2fc = direct_activated$log2FoldChange,
            stringsAsFactors = FALSE
        )
        edges_rep <- data.frame(
            from = chip_base,
            to = as.character(direct_repressed[[peak_gene_col]]),
            type = "repression",
            weight = if (!is.null(score_col)) direct_repressed[[score_col]] else 1,
            log2fc = direct_repressed$log2FoldChange,
            stringsAsFactors = FALSE
        )
        edges_all <- rbind(edges_act, edges_rep)

        if (nrow(edges_all) > 0L) {
            g <- igraph::graph_from_data_frame(edges_all, directed = TRUE)
            igraph::V(g)$type <- ifelse(igraph::V(g)$name == chip_base, "TF", "target")
            igraph::V(g)$color <- ifelse(igraph::V(g)$name == chip_base, "red",
                ifelse(igraph::V(g)$name %in% edges_act$to, "green", "blue")
            )

            graphml_path <- file.path(net_dir, sprintf("%s_network.graphml", run_label))
            igraph::write_graph(g, graphml_path, format = "graphml")
            message(sprintf(
                "[network] %d nodes, %d edges → %s",
                igraph::vcount(g), igraph::ecount(g), graphml_path
            ))

            deg_centrality <- igraph::degree(g, mode = "in")
            top_hubs <- sort(deg_centrality[names(deg_centrality) != chip_base],
                decreasing = TRUE
            )[seq_len(min(20L, length(deg_centrality) - 1L))]

            write_tsv(
                data.frame(
                    gene_id = names(top_hubs), in_degree = as.integer(top_hubs),
                    stringsAsFactors = FALSE
                ),
                file.path(out_dir, "network_hub_genes.tsv")
            )
        }
    },
    error = function(e) message("[warn] Network construction skipped: ", conditionMessage(e))
)

message("[validation] Checking for known p53 targets")

known_p53_targets <- c(
    "CDKN1A", "BAX", "MDM2", "BBC3", "GADD45A",
    "TIGAR", "GDF15", "PUMA", "TP53I3", "SESN1"
)

sym_col <- if ("SYMBOL" %in% colnames(direct_targets)) "SYMBOL" else if ("geneSymbol" %in% colnames(direct_targets)) "geneSymbol" else NULL

if (!is.null(sym_col)) {
    validated <- direct_targets[
        !is.na(direct_targets[[sym_col]]) & direct_targets[[sym_col]] %in% known_p53_targets,
        c(
            peak_gene_col, sym_col, "log2FoldChange", "padj",
            if (!is.null(score_col)) score_col else character(0)
        )
    ]
    message(sprintf(
        "[validation] Known p53 targets found: %d / %d",
        nrow(validated), length(known_p53_targets)
    ))
    if (nrow(validated) > 0L) {
        message(paste(" ", validated[[sym_col]], collapse = "\n"))
        write_tsv(validated, file.path(out_dir, "known_p53_targets_validated.tsv"))
    }
} else {
    message("[validation] Symbol column not found — skipping known target check")
}

message(sprintf("[done] Integration tables → %s", out_dir))
message(sprintf("[done] Network           → %s", net_dir))
message(sprintf("[done] Figures           → %s", fig_dir))
