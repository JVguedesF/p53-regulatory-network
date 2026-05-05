#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(tximport)
    library(DESeq2)
    library(ggplot2)
    library(org.Hs.eg.db)
    library(AnnotationDbi)
    library(ComplexHeatmap)
    library(grid)
})

FASTA_GZ <- "data/genome/fasta/Homo_sapiens.GRCh38.cdna.all.fa.gz"
TX2GENE_CACHE <- "data/genome/tx2gene.tsv"
MIN_COUNT <- 10L
MIN_SAMPLE_FRAC <- 0.25

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
    stop(
        "Usage: Rscript 09_differential_expression.R <metadata.tsv>",
        " [treated_condition] [reference_condition]",
        " [lfc_threshold=1] [padj_threshold=0.05]",
        call. = FALSE
    )
}

tsv_path <- args[1]
treated_cond <- if (length(args) >= 2) args[2] else NULL
reference_cond <- if (length(args) >= 3) args[3] else NULL
lfc_threshold <- if (length(args) >= 4) as.numeric(args[4]) else 1.0
padj_threshold <- if (length(args) >= 5) as.numeric(args[5]) else 0.05

tsv_base <- sub("\\.tsv$", "", basename(tsv_path))

save_pdf <- function(expr, path, width = 8, height = 6) {
    pdf(path, width = width, height = height)
    on.exit(dev.off())
    force(expr)
}

mark_done <- function(tsv, accession) {
    system2(
        "python",
        args   = c("src/tsv_updater.py", tsv, accession, "DEG_Status=DONE"),
        stdout = FALSE,
        stderr = FALSE
    )
}

build_tx2gene <- function(fasta_gz, cache_path) {
    if (file.exists(cache_path)) {
        message("[tx2gene] Loading cached mapping")
        df <- read.table(
            cache_path,
            header           = FALSE,
            sep              = "\t",
            col.names        = c("tx_id", "gene_id"),
            stringsAsFactors = FALSE
        )
        message(sprintf(
            "[tx2gene] %d transcripts → %d genes",
            nrow(df), length(unique(df$gene_id))
        ))
        return(df)
    }

    if (!file.exists(fasta_gz)) {
        stop("cDNA FASTA not found: ", fasta_gz,
            "\nExpected from script 08 Salmon index build.",
            call. = FALSE
        )
    }

    message("[tx2gene] Building from cDNA FASTA (cached after this run)...")
    dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)

    awk_prog <- paste0(
        "/^>/{",
        "tx=$1; sub(/^>/,\"\",tx); sub(/\\.[0-9]+$/,\"\",tx);",
        "g=\"\";",
        "for(i=2;i<=NF;i++) if($i~/^gene:/){g=$i;",
        "  sub(/^gene:/,\"\",g); sub(/\\.[0-9]+$/,\"\",g)};",
        "if(g!=\"\") print tx\"\\t\"g}"
    )

    ret <- system2("bash", args = c("-c", paste(
        "gzip -dc", shQuote(fasta_gz),
        "| awk",    shQuote(awk_prog),
        ">",        shQuote(cache_path)
    )))

    if (ret != 0 || !file.exists(cache_path) || file.size(cache_path) == 0) {
        stop("Failed to build tx2gene from FASTA. Check file integrity.", call. = FALSE)
    }

    df <- read.table(
        cache_path,
        header           = FALSE,
        sep              = "\t",
        col.names        = c("tx_id", "gene_id"),
        stringsAsFactors = FALSE
    )
    message(sprintf(
        "[tx2gene] %d transcripts → %d genes",
        nrow(df), length(unique(df$gene_id))
    ))
    df
}

if (!file.exists(tsv_path)) stop("TSV not found: ", tsv_path, call. = FALSE)

meta <- read.table(
    tsv_path,
    sep              = "\t",
    header           = TRUE,
    stringsAsFactors = FALSE,
    quote            = "",
    fill             = TRUE,
    comment.char     = ""
)

rna_meta <- meta[
    toupper(meta$Sample_Type) == "RNA" &
        !is.na(meta$BAM_Path) &
        nzchar(meta$BAM_Path),
]

if (nrow(rna_meta) == 0L) {
    stop("No RNA rows with quantification paths found.\nRun script 08 first.", call. = FALSE)
}

samples <- rna_meta[!duplicated(rna_meta$BAM_Path), ]

missing_sf <- !sapply(samples$BAM_Path, function(d) file.exists(file.path(d, "quant.sf")))
if (any(missing_sf)) {
    stop("Missing quant.sf in:\n",
        paste(samples$BAM_Path[missing_sf], collapse = "\n"),
        "\nRun script 08 first.",
        call. = FALSE
    )
}

conditions <- unique(samples$Condition)

if (is.null(treated_cond)) {
    if (length(conditions) != 2L) {
        stop(sprintf(
            "Cannot auto-detect conditions: %d found (%s).\nSpecify: Rscript 09_... <tsv> <treated> <reference>",
            length(conditions), paste(conditions, collapse = ", ")
        ), call. = FALSE)
    }
    treated_cond <- conditions[1]
    reference_cond <- conditions[2]
    message(sprintf("[auto] treated=%s  reference=%s", treated_cond, reference_cond))
}

for (cond in c(treated_cond, reference_cond)) {
    if (!cond %in% conditions) {
        stop(sprintf(
            "Condition '%s' not found. Available: %s",
            cond, paste(conditions, collapse = ", ")
        ), call. = FALSE)
    }
}

samples <- samples[samples$Condition %in% c(treated_cond, reference_cond), ]
samples$condition <- factor(samples$Condition, levels = c(reference_cond, treated_cond))
rownames(samples) <- samples$Accession

message(sprintf(
    "[setup] %s | %s vs %s | %d samples (%d treated, %d reference)",
    tsv_base, treated_cond, reference_cond, nrow(samples),
    sum(samples$Condition == treated_cond),
    sum(samples$Condition == reference_cond)
))

if (any(table(samples$condition) < 2L)) {
    warning("Some conditions have fewer than 2 replicates. Dispersion estimation will be unreliable.")
}

tx2gene <- build_tx2gene(FASTA_GZ, TX2GENE_CACHE)

quant_files <- file.path(samples$BAM_Path, "quant.sf")
names(quant_files) <- samples$Accession

message(sprintf("[tximport] Importing %d sample(s)", length(quant_files)))
txi <- tximport(
    quant_files,
    type            = "salmon",
    tx2gene         = tx2gene,
    ignoreTxVersion = TRUE
)

dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~condition)

min_samples <- ceiling(ncol(dds) * MIN_SAMPLE_FRAC)
dds <- dds[rowSums(counts(dds) >= MIN_COUNT) >= min_samples, ]

message(sprintf("[DESeq2] Running on %d genes (after low-count filter)", nrow(dds)))
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", treated_cond, reference_cond), alpha = padj_threshold)
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

message("[DESeq2] Summary:")
print(summary(res))

sym_map <- tryCatch(
    AnnotationDbi::select(
        org.Hs.eg.db,
        keys    = rownames(res_df),
        columns = c("SYMBOL", "GENENAME"),
        keytype = "ENSEMBL"
    ),
    error = function(e) {
        message("[warn] Symbol annotation failed: ", conditionMessage(e))
        data.frame(ENSEMBL = rownames(res_df), SYMBOL = NA_character_, GENENAME = NA_character_)
    }
)
sym_map <- sym_map[!duplicated(sym_map$ENSEMBL), ]
res_df <- merge(res_df, sym_map, by.x = "gene_id", by.y = "ENSEMBL", all.x = TRUE)
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]

label <- sprintf("%s_vs_%s", treated_cond, reference_cond)

# DESeq2 outputs feed script 10 → pipeline_outputs
out_dir <- file.path("pipeline_outputs", "rnaseq", tsv_base, "deseq2", label)
deg_dir <- file.path("pipeline_outputs", "rnaseq", tsv_base, "deseq2")

# gene lists and figures are for humans → results/
tab_dir <- file.path("results", "tables", "rnaseq", tsv_base)
fig_dir <- file.path("results", "figures", "rnaseq", tsv_base)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

write.table(res_df,
    file.path(out_dir, "degs_all.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)

degs_sig <- res_df[
    !is.na(res_df$padj) &
        res_df$padj < padj_threshold &
        abs(res_df$log2FoldChange) > lfc_threshold,
]

write.table(degs_sig,
    file.path(out_dir, "degs_significant.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)

n_up <- sum(degs_sig$log2FoldChange > 0, na.rm = TRUE)
n_down <- sum(degs_sig$log2FoldChange < 0, na.rm = TRUE)
message(sprintf("[DESeq2] DEGs: %d total | %d up | %d down", nrow(degs_sig), n_up, n_down))

activated <- degs_sig$gene_id[!is.na(degs_sig$log2FoldChange) & degs_sig$log2FoldChange > 0]
repressed <- degs_sig$gene_id[!is.na(degs_sig$log2FoldChange) & degs_sig$log2FoldChange < 0]

writeLines(activated, file.path(tab_dir, "activated_genes.txt"))
writeLines(repressed, file.path(tab_dir, "repressed_genes.txt"))

vsd <- vst(dds, blind = TRUE)

save_pdf(
    plotMA(res,
        main  = sprintf("MA Plot — %s vs %s", treated_cond, reference_cond),
        ylim  = c(-5, 5),
        alpha = padj_threshold
    ),
    file.path(fig_dir, "ma_plot.pdf")
)

vol_df <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]
vol_df$status <- "NS"
vol_df$status[vol_df$padj < padj_threshold & vol_df$log2FoldChange > lfc_threshold] <- "Up"
vol_df$status[vol_df$padj < padj_threshold & vol_df$log2FoldChange < -lfc_threshold] <- "Down"
vol_df$status <- factor(vol_df$status, levels = c("Up", "Down", "NS"))

save_pdf(
    print(
        ggplot(vol_df, aes(x = log2FoldChange, y = -log10(padj), color = status)) +
            geom_point(alpha = 0.5, size = 0.8, shape = 16) +
            scale_color_manual(
                values = c(Up = "#E74C3C", Down = "#2980B9", NS = "grey70"),
                drop   = FALSE
            ) +
            geom_vline(
                xintercept = c(-lfc_threshold, lfc_threshold),
                linetype = "dashed", linewidth = 0.4, color = "black"
            ) +
            geom_hline(
                yintercept = -log10(padj_threshold),
                linetype = "dashed", linewidth = 0.4, color = "black"
            ) +
            labs(
                title = sprintf("Volcano — %s vs %s", treated_cond, reference_cond),
                x     = "Log2 Fold Change",
                y     = "-Log10(adjusted p-value)",
                color = NULL
            ) +
            theme_bw(base_size = 12) +
            theme(legend.position = "top")
    ),
    file.path(fig_dir, "volcano_plot.pdf")
)

save_pdf(
    {
        pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
        pct_var <- round(100 * attr(pca_data, "percentVar"))
        print(
            ggplot(pca_data, aes(PC1, PC2, color = condition)) +
                geom_point(size = 3) +
                geom_text(aes(label = name), size = 2.5, vjust = -0.9, show.legend = FALSE) +
                xlab(sprintf("PC1: %d%% variance", pct_var[1])) +
                ylab(sprintf("PC2: %d%% variance", pct_var[2])) +
                scale_color_manual(
                    values = setNames(c("#2980B9", "#E74C3C"), c(reference_cond, treated_cond))
                ) +
                ggtitle("PCA — Variance-Stabilized Counts") +
                theme_bw(base_size = 12)
        )
    },
    file.path(fig_dir, "pca_plot.pdf")
)

if (nrow(degs_sig) >= 2L) {
    top_n <- min(50L, nrow(degs_sig))
    top_ids <- degs_sig$gene_id[seq_len(top_n)]
    top_ids_present <- intersect(top_ids, rownames(assay(vsd)))

    mat <- assay(vsd)[top_ids_present, , drop = FALSE]
    mat <- t(scale(t(mat)))

    sym_lookup <- setNames(
        degs_sig$SYMBOL[match(rownames(mat), degs_sig$gene_id)],
        rownames(mat)
    )
    row_labels <- ifelse(!is.na(sym_lookup) & nzchar(sym_lookup), sym_lookup, rownames(mat))

    col_ann <- HeatmapAnnotation(
        condition = as.character(samples[colnames(mat), "condition"]),
        col = list(condition = setNames(
            c("#2980B9", "#E74C3C"),
            c(reference_cond, treated_cond)
        ))
    )

    save_pdf(
        draw(Heatmap(
            mat,
            name = "Z-score",
            top_annotation = col_ann,
            row_labels = row_labels,
            row_names_gp = gpar(fontsize = 7),
            cluster_rows = TRUE,
            cluster_columns = FALSE,
            show_column_names = TRUE,
            column_title = sprintf(
                "Top %d DEGs — %s vs %s", nrow(mat), treated_cond, reference_cond
            )
        )),
        file.path(fig_dir, "heatmap_top50.pdf"),
        width = 8, height = 10
    )
}

for (acc in rna_meta$Accession) mark_done(tsv_path, acc)

message(sprintf("[done] DESeq2 tables → %s", out_dir))
message(sprintf("[done] Gene lists   → %s", tab_dir))
message(sprintf("[done] Figures      → %s", fig_dir))
