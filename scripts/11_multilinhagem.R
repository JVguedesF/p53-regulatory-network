suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(ggupset)
    library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)

integration_root <- if (length(args) >= 1) args[1] else file.path("pipeline_outputs", "integration")
lfc_threshold <- if (length(args) >= 2) as.numeric(args[2]) else 1.0
padj_threshold <- if (length(args) >= 3) as.numeric(args[3]) else 0.05

out_dir <- file.path("pipeline_outputs", "multilinhagem")
fig_dir <- file.path("results", "figures", "multilinhagem")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv_safe <- function(path) {
    read.table(path,
        sep = "\t", header = TRUE, stringsAsFactors = FALSE,
        quote = "", fill = TRUE, comment.char = ""
    )
}

save_pdf <- function(expr, path, width = 8, height = 6) {
    pdf(NULL)
    pdf(path, width = width, height = height)
    on.exit(dev.off())
    force(expr)
}

get_sym_col <- function(df) {
    candidates <- c("SYMBOL", "geneSymbol")
    found <- intersect(candidates, colnames(df))
    if (length(found) == 0L) NULL else found[1]
}

target_files <- list.files(
    integration_root,
    pattern      = "direct_targets_all\\.tsv$",
    recursive    = TRUE,
    full.names   = TRUE
)

if (length(target_files) == 0L) {
    stop(
        "No direct_targets_all.tsv files found under: ", integration_root,
        "\nRun 09_integration.R for each ChIP+RNA pair first.",
        call. = FALSE
    )
}

message(sprintf("[setup] Found %d integration run(s)", length(target_files)))

runs <- setNames(target_files, basename(dirname(target_files)))

targets_long <- do.call(rbind, lapply(names(runs), function(run_label) {
    df <- read_tsv_safe(runs[[run_label]])
    sym_col <- get_sym_col(df)

    gene_id_col <- if ("gene_id" %in% colnames(df)) "gene_id" else colnames(df)[1]
    gene_sym <- if (!is.null(sym_col)) df[[sym_col]] else df[[gene_id_col]]

    data.frame(
        run = run_label,
        gene_id = df[[gene_id_col]],
        symbol = gene_sym,
        log2FC = df$log2FoldChange,
        padj = df$padj,
        direction = ifelse(df$log2FoldChange > 0, "activated", "repressed"),
        stringsAsFactors = FALSE
    )
}))

message(sprintf("[load] Total direct-target records across all runs: %d", nrow(targets_long)))

all_symbols <- sort(unique(targets_long$symbol[!is.na(targets_long$symbol)]))
run_labels <- names(runs)

conservation_mat <- matrix(
    NA_character_,
    nrow = length(all_symbols),
    ncol = length(run_labels),
    dimnames = list(all_symbols, run_labels)
)

for (i in seq_len(nrow(targets_long))) {
    sym <- targets_long$symbol[i]
    run <- targets_long$run[i]
    if (!is.na(sym) && sym %in% all_symbols && run %in% run_labels) {
        conservation_mat[sym, run] <- targets_long$direction[i]
    }
}

n_runs_per_gene <- rowSums(!is.na(conservation_mat))

conserved_targets <- names(n_runs_per_gene[n_runs_per_gene == length(run_labels)])
message(sprintf(
    "[conservation] Targets shared across ALL %d run(s): %d",
    length(run_labels), length(conserved_targets)
))

context_specific <- names(n_runs_per_gene[n_runs_per_gene == 1L])
message(sprintf("[conservation] Context-specific targets (1 run only): %d", length(context_specific)))

known_p53 <- c(
    "CDKN1A", "BAX", "MDM2", "BBC3", "GADD45A",
    "TIGAR", "GDF15", "PUMA", "TP53I3", "SESN1"
)

known_found <- intersect(known_p53, all_symbols)
known_conserved <- intersect(known_p53, conserved_targets)

message(sprintf(
    "[validation] Known p53 targets detected: %d / %d  |  conserved across all runs: %d",
    length(known_found), length(known_p53), length(known_conserved)
))

conservation_df <- data.frame(
    symbol = rownames(conservation_mat),
    n_runs = n_runs_per_gene,
    known_p53 = rownames(conservation_mat) %in% known_p53,
    conservation_mat,
    check.names = FALSE,
    stringsAsFactors = FALSE
)
conservation_df <- conservation_df[order(-conservation_df$n_runs, conservation_df$symbol), ]

write.table(conservation_df,
    file.path(out_dir, "target_conservation_matrix.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
    conservation_df[conservation_df$symbol %in% conserved_targets, ],
    file.path(out_dir, "conserved_targets_all_runs.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
    conservation_df[conservation_df$symbol %in% context_specific, ],
    file.path(out_dir, "context_specific_targets.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)

run_summary <- do.call(rbind, lapply(run_labels, function(r) {
    sub_df <- targets_long[targets_long$run == r, ]
    data.frame(
        run = r,
        n_total = nrow(sub_df),
        n_act = sum(sub_df$direction == "activated", na.rm = TRUE),
        n_rep = sum(sub_df$direction == "repressed", na.rm = TRUE),
        n_known_p53 = sum(sub_df$symbol %in% known_p53, na.rm = TRUE),
        stringsAsFactors = FALSE
    )
}))

write.table(run_summary,
    file.path(out_dir, "per_run_summary.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)

message("[plots] Bar chart: targets per run")

bar_df <- targets_long %>%
    count(run, direction) %>%
    mutate(run = factor(run, levels = run_labels))

save_pdf(
    print(
        ggplot(bar_df, aes(x = run, y = n, fill = direction)) +
            geom_col(position = "stack", width = 0.6) +
            scale_fill_manual(values = c(activated = "#E74C3C", repressed = "#2980B9")) +
            labs(
                title = "Direct p53 Targets per Integration Run",
                x     = NULL,
                y     = "Number of genes",
                fill  = NULL
            ) +
            theme_bw(base_size = 12) +
            theme(
                axis.text.x = element_text(angle = 30, hjust = 1),
                legend.position = "top"
            )
    ),
    file.path(fig_dir, "targets_per_run.pdf"),
    width = max(5, 2 * length(run_labels)), height = 5
)

message("[plots] Histogram: conservation across runs")

hist_df <- data.frame(n_runs = n_runs_per_gene)

save_pdf(
    print(
        ggplot(hist_df, aes(x = n_runs)) +
            geom_histogram(
                binwidth = 1, fill = "#2ECC71", color = "white", boundary = 0.5
            ) +
            scale_x_continuous(breaks = seq_len(length(run_labels))) +
            labs(
                title = "Target Gene Conservation across Runs",
                x     = "Number of runs gene is a direct target",
                y     = "Gene count"
            ) +
            theme_bw(base_size = 12)
    ),
    file.path(fig_dir, "conservation_histogram.pdf"),
    width = 6, height = 4
)

if (length(run_labels) >= 2L) {
    message("[plots] UpSet plot")

    upset_df <- do.call(rbind, lapply(all_symbols, function(sym) {
        present_runs <- run_labels[!is.na(conservation_mat[sym, ])]
        data.frame(symbol = sym, runs = I(list(present_runs)), stringsAsFactors = FALSE)
    }))

    save_pdf(
        print(
            ggplot(upset_df, aes(x = runs)) +
                geom_bar(fill = "#8E44AD") +
                scale_x_upset(order_by = "freq") +
                labs(
                    title = "UpSet: Overlap of Direct p53 Targets",
                    x     = NULL,
                    y     = "Gene count"
                ) +
                theme_bw(base_size = 12) +
                theme_combmatrix(combmatrix.label.text = element_text(size = 9))
        ),
        file.path(fig_dir, "upset_targets.pdf"),
        width = max(6, 2 * length(run_labels)), height = 5
    )
}

if (length(run_labels) >= 2L && length(conserved_targets) > 0L) {
    message("[plots] Heatmap: log2FC of conserved targets")

    lfc_mat <- matrix(
        NA_real_,
        nrow = length(conserved_targets),
        ncol = length(run_labels),
        dimnames = list(conserved_targets, run_labels)
    )

    for (run in run_labels) {
        sub <- targets_long[targets_long$run == run & targets_long$symbol %in% conserved_targets, ]
        lfc_mat[sub$symbol, run] <- sub$log2FC
    }

    max_rows <- 40L
    mean_abs <- rowMeans(abs(lfc_mat), na.rm = TRUE)
    keep_genes <- names(sort(mean_abs, decreasing = TRUE))[seq_len(min(max_rows, nrow(lfc_mat)))]
    lfc_mat <- lfc_mat[keep_genes, , drop = FALSE]

    row_anno <- data.frame(
        Known_p53 = ifelse(rownames(lfc_mat) %in% known_p53, "yes", "no"),
        row.names = rownames(lfc_mat)
    )

    save_pdf(
        pheatmap(
            lfc_mat,
            annotation_row = row_anno,
            annotation_colors = list(Known_p53 = c(yes = "#F39C12", no = "grey85")),
            color = colorRampPalette(c("#2980B9", "white", "#E74C3C"))(101),
            breaks = seq(-3, 3, length.out = 102),
            cluster_cols = length(run_labels) > 1L,
            cluster_rows = TRUE,
            border_color = NA,
            fontsize_row = 8,
            main = sprintf(
                "log2FC of Conserved p53 Targets (top %d by mean |lfc|)", max_rows
            ),
            silent = TRUE
        ),
        file.path(fig_dir, "heatmap_conserved_lfc.pdf"),
        width = max(5, 2 * length(run_labels) + 2),
        height = max(6, 0.18 * nrow(lfc_mat) + 3)
    )
}

message(sprintf("[done] Tables  → %s", out_dir))
message(sprintf("[done] Figures → %s", fig_dir))
