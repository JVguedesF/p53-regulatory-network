#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
    stop(
        "Usage: Rscript 06_peak_annotation.R <metadata.tsv>",
        " [tss_upstream] [tss_downstream]",
        call. = FALSE
    )
}

tsv_path <- args[1]
tss_upstream <- if (length(args) >= 2) as.integer(args[2]) else 3000L
tss_downstream <- if (length(args) >= 3) as.integer(args[3]) else 3000L

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

save_pdf <- function(expr, path, width = 7, height = 6) {
    pdf(path, width = width, height = height)
    on.exit(dev.off())
    force(expr)
}

mark_annotated <- function(tsv, accession) {
    system2(
        "python",
        args   = c("src/tsv_updater.py", tsv, accession, "Annotated=TRUE"),
        stdout = FALSE,
        stderr = FALSE
    )
}

if (!file.exists(tsv_path)) {
    stop("TSV not found: ", tsv_path, call. = FALSE)
}

meta <- read.table(
    tsv_path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE
)

samples <- meta[
    meta$Sample_Type == "IP" &
        !is.na(meta$Peak_File) &
        nzchar(meta$Peak_File),
]

if (nrow(samples) == 0) {
    stop("No IP rows with Peak_File found in TSV.", call. = FALSE)
}

message(sprintf("[setup] %d IP sample(s) to process", nrow(samples)))

anno_list <- list()

for (i in seq_len(nrow(samples))) {
    row <- samples[i, ]
    acc <- row$Accession
    label <- sprintf("%s_Rep%s_%s", row$Condition, row$Replicate, acc)

    sample_dir <- file.path("results", "annotation", label)
    rdata_path <- file.path(sample_dir, paste0(label, "_anno.RData"))

    if (isTRUE(row$Annotated == "TRUE") && file.exists(rdata_path)) {
        message(sprintf("[skip] %s — already annotated", label))
        load(rdata_path)
        anno_list[[label]] <- peak_anno
        next
    }

    if (!file.exists(row$Peak_File)) {
        message(sprintf("[warn] Peak file not found for %s — skipping", label))
        next
    }

    message(sprintf("[annotate] %s", label))
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)

    peaks <- readPeakFile(row$Peak_File, as = "GRanges")
    peak_anno <- annotatePeak(
        peaks,
        tssRegion = c(-tss_upstream, tss_downstream),
        TxDb      = txdb,
        annoDb    = "org.Hs.eg.db",
        verbose   = FALSE
    )

    anno_df <- as.data.frame(peak_anno)
    write.table(
        anno_df,
        file.path(sample_dir, paste0(label, "_annotation.tsv")),
        sep       = "\t",
        quote     = FALSE,
        row.names = FALSE
    )

    save_pdf(
        plotAnnoPie(peak_anno),
        file.path(sample_dir, paste0(label, "_annopie.pdf"))
    )
    save_pdf(
        plotAnnoBar(peak_anno),
        file.path(sample_dir, paste0(label, "_annobar.pdf"))
    )
    save_pdf(
        plotDistToTSS(peak_anno),
        file.path(sample_dir, paste0(label, "_dist_to_tss.pdf"))
    )

    save(peak_anno, file = rdata_path)
    mark_annotated(tsv_path, acc)

    anno_list[[label]] <- peak_anno
    message(sprintf("[done] %s — %d peaks", label, length(peaks)))
}

if (length(anno_list) == 0) {
    stop("No samples were successfully annotated.", call. = FALSE)
}

message("[comparative] Generating multi-sample AnnoBar")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

save_pdf(
    plotAnnoBar(anno_list),
    "results/figures/annobar_comparative.pdf",
    width  = 9,
    height = 5
)

message("[avgprof] Computing tag matrices")

promoter_windows <- getPromoters(
    TxDb       = txdb,
    upstream   = tss_upstream,
    downstream = tss_downstream
)

tag_matrix_list <- lapply(anno_list, function(anno) {
    getTagMatrix(as(anno, "GRanges"), windows = promoter_windows)
})

save_pdf(
    plotAvgProf(
        tag_matrix_list,
        xlim = c(-tss_upstream, tss_downstream),
        conf = 0.95,
        resample = 500,
        xlab = "Distance to TSS (bp)",
        ylab = "Read count frequency"
    ),
    "results/figures/avgprof_tss.pdf",
    width = 8,
    height = 5
)

message("[multi-peak] Identifying genes with multiple promoter peaks")

all_anno_df <- do.call(rbind, lapply(names(anno_list), function(lbl) {
    df <- as.data.frame(anno_list[[lbl]])
    df$label <- lbl
    df
}))

promoter_df <- all_anno_df[grepl("Promoter", all_anno_df$annotation), ]
peak_counts <- sort(table(promoter_df$geneId), decreasing = TRUE)
multi_peak <- peak_counts[peak_counts > 1]

if (length(multi_peak) > 0) {
    message(sprintf(
        "[multi-peak] %d gene(s) with >1 promoter peak",
        length(multi_peak)
    ))
    print(head(multi_peak, 20))
} else {
    message("[multi-peak] No genes with multiple promoter peaks found")
}

dir.create("results/peaks", recursive = TRUE, showWarnings = FALSE)

promoter_genes <- unique(promoter_df$geneId)
writeLines(promoter_genes, "results/peaks/p53_promoter_genes.txt")

message(sprintf(
    "[export] %d unique promoter genes → results/peaks/p53_promoter_genes.txt",
    length(promoter_genes)
))
