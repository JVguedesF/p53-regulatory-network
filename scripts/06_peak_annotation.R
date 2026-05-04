#!/usr/bin/env Rscript

pdf(NULL)

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
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
max_peaks_prof <- 5000L

tsv_base <- sub("\\.tsv$", "", basename(tsv_path))
main_chrs <- paste0("chr", c(1:22, "X", "Y"))

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
GenomeInfoDb::seqlevels(txdb, pruning.mode = "coarse") <- main_chrs

promoter_windows <- getPromoters(
    TxDb       = txdb,
    upstream   = tss_upstream,
    downstream = tss_downstream
)

save_pdf <- function(expr, path, width = 7, height = 6) {
    pdf(path, width = width, height = height)
    on.exit(dev.off())
    force(expr)
}

read_tsv <- function(path) {
    read.table(
        path,
        sep              = "\t",
        header           = TRUE,
        stringsAsFactors = FALSE,
        quote            = "",
        fill             = TRUE,
        comment.char     = ""
    )
}

mark_annotated <- function(tsv, accession) {
    system2(
        "python",
        args   = c("src/tsv_updater.py", tsv, accession, "Annotated=TRUE"),
        stdout = FALSE,
        stderr = FALSE
    )
}

process_sample <- function(row, tsv_path, tsv_base, txdb, promoter_windows,
                           main_chrs, max_peaks_prof,
                           tss_upstream, tss_downstream) {
    acc <- row$Accession
    label <- sprintf("%s_Rep%s_%s", row$Condition, row$Replicate, acc)
    sample_dir <- file.path("results", "annotation", tsv_base, label)
    rdata_path <- file.path(sample_dir, paste0(label, "_anno.RData"))

    if (isTRUE(row$Annotated == "TRUE") && file.exists(rdata_path)) {
        message(sprintf("[skip] %s — already annotated", label))
        return(label)
    }

    if (!file.exists(row$Peak_File)) {
        message(sprintf("[warn] Peak file not found for %s — skipping", label))
        return(NULL)
    }

    message(sprintf("[annotate] %s", label))
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)

    peaks <- readPeakFile(row$Peak_File, as = "GRanges")

    if (!any(grepl("^chr", GenomeInfoDb::seqlevels(peaks)))) {
        GenomeInfoDb::seqlevels(peaks) <- paste0("chr", GenomeInfoDb::seqlevels(peaks))
    }

    peaks <- GenomeInfoDb::keepSeqlevels(
        peaks,
        intersect(GenomeInfoDb::seqlevels(peaks), main_chrs),
        pruning.mode = "coarse"
    )

    peak_anno <- suppressWarnings(
        annotatePeak(
            peaks,
            tssRegion = c(-tss_upstream, tss_downstream),
            TxDb      = txdb,
            annoDb    = "org.Hs.eg.db",
            verbose   = FALSE
        )
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

    rm(anno_df)
    gc()

    n_peaks <- length(peaks)
    gr <- peaks

    if (length(gr) > max_peaks_prof) {
        score_col <- if ("V5" %in% colnames(mcols(gr))) "V5" else colnames(mcols(gr))[1]
        top_idx <- order(mcols(gr)[[score_col]], decreasing = TRUE)[seq_len(max_peaks_prof)]
        gr <- gr[top_idx]
        message(sprintf("[avgprof] %s — downsampled to top %d peaks by score", label, max_peaks_prof))
    }

    if (length(gr) == 0) {
        message(sprintf("[avgprof] %s — no peaks on main chromosomes, skipping", label))
    } else {
        tag_matrix <- tryCatch(
            getTagMatrix(gr, windows = promoter_windows),
            error = function(e) {
                message(sprintf("[avgprof] %s — getTagMatrix failed: %s", label, conditionMessage(e)))
                NULL
            }
        )

        if (!is.null(tag_matrix)) {
            tag_matrix <- tag_matrix[rowSums(tag_matrix) > 0, , drop = FALSE]

            if (nrow(tag_matrix) < 2) {
                message(sprintf("[avgprof] %s — fewer than 2 non-zero rows, skipping avgprof", label))
            } else {
                tryCatch(
                    save_pdf(
                        plotAvgProf(
                            tag_matrix,
                            xlim  = c(-tss_upstream, tss_downstream),
                            conf  = FALSE,
                            xlab  = "Distance to TSS (bp)",
                            ylab  = "Read count frequency"
                        ),
                        file.path(sample_dir, paste0(label, "_avgprof.pdf")),
                        width = 8,
                        height = 5
                    ),
                    error = function(e) {
                        message(sprintf("[avgprof] %s — plotAvgProf failed: %s", label, conditionMessage(e)))
                    }
                )
            }

            rm(tag_matrix)
            gc()
        }
    }

    rm(gr, peaks)
    gc()

    save(peak_anno, file = rdata_path)
    mark_annotated(tsv_path, acc)

    message(sprintf("[done] %s — %d peaks", label, n_peaks))

    rm(peak_anno)
    gc()

    return(label)
}

if (!file.exists(tsv_path)) {
    stop("TSV not found: ", tsv_path, call. = FALSE)
}

meta <- read_tsv(tsv_path)
samples <- meta[
    meta$Sample_Type == "IP" &
        !is.na(meta$Peak_File) &
        nzchar(meta$Peak_File),
]

if (nrow(samples) == 0) {
    stop("No IP rows with Peak_File found in TSV.", call. = FALSE)
}

message(sprintf("[setup] %s | %d IP sample(s) to process", tsv_base, nrow(samples)))

processed_labels <- c()

for (i in seq_len(nrow(samples))) {
    label <- process_sample(
        row              = samples[i, ],
        tsv_path         = tsv_path,
        tsv_base         = tsv_base,
        txdb             = txdb,
        promoter_windows = promoter_windows,
        main_chrs        = main_chrs,
        max_peaks_prof   = max_peaks_prof,
        tss_upstream     = tss_upstream,
        tss_downstream   = tss_downstream
    )

    if (!is.null(label)) {
        processed_labels <- c(processed_labels, label)
    }

    gc()
}

if (length(processed_labels) == 0) {
    stop("No samples were successfully annotated.", call. = FALSE)
}

dataset_dir <- file.path("results", "annotation", tsv_base)
dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)

message("[comparative] Generating multi-sample AnnoBar")

anno_list <- list()

for (label in processed_labels) {
    rdata_path <- file.path(dataset_dir, label, paste0(label, "_anno.RData"))
    if (file.exists(rdata_path)) {
        local({
            load(rdata_path)
            anno_list[[label]] <<- peak_anno
        })
    }
    gc()
}

if (length(anno_list) > 1) {
    save_pdf(
        plotAnnoBar(anno_list),
        file.path(dataset_dir, "annobar_comparative.pdf"),
        width  = 9,
        height = 5
    )
}

rm(anno_list)
gc()

message("[multi-peak] Identifying genes with multiple promoter peaks")

all_anno_df <- do.call(rbind, lapply(processed_labels, function(label) {
    tsv_file <- file.path(dataset_dir, label, paste0(label, "_annotation.tsv"))
    if (!file.exists(tsv_file)) {
        return(NULL)
    }
    df <- read_tsv(tsv_file)
    df$label <- label
    df
}))

promoter_df <- all_anno_df[grepl("Promoter", all_anno_df$annotation), ]
peak_counts <- sort(table(promoter_df$geneId), decreasing = TRUE)
multi_peak <- peak_counts[peak_counts > 1]

if (length(multi_peak) > 0) {
    message(sprintf("[multi-peak] %d gene(s) with >1 promoter peak", length(multi_peak)))
    print(head(multi_peak, 20))
} else {
    message("[multi-peak] No genes with multiple promoter peaks found")
}

rm(all_anno_df)
gc()

promoter_genes <- as.character(unique(promoter_df$geneId))
writeLines(promoter_genes, file.path(dataset_dir, "promoter_genes.txt"))

message(sprintf(
    "[export] %d unique promoter genes → %s/promoter_genes.txt",
    length(promoter_genes), dataset_dir
))
