#!/usr/bin/env Rscript

pdf(NULL)

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(GenomeInfoDb)
    library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
    stop(
        "Usage: Rscript scripts/07_peak_annotation.R <metadata.tsv> [<metadata2.tsv> ...]",
        call. = FALSE
    )
}

tss_upstream   <- 3000L
tss_downstream <- 3000L
max_peaks_prof <- 2000L

tsv_paths <- args

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
main_chrs <- paste0("chr", c(1:22, "X", "Y"))
GenomeInfoDb::seqlevels(txdb, pruning.mode = "coarse") <- main_chrs

promoter_windows <- getPromoters(
    TxDb       = txdb,
    upstream   = tss_upstream,
    downstream = tss_downstream
)

sanitize_id <- function(x) {
    x <- gsub("[[:space:]]+", "_", x)
    x <- gsub("[^A-Za-z0-9_.+-]+", "_", x)
    gsub("^_+|_+$", "", x)
}

sample_id <- function(condition, sample_type, replicate) {
    sanitize_id(sprintf("%s_%s_Rep%s", condition, sample_type, replicate))
}

save_pdf <- function(expr, path, width = 7, height = 6) {
    pdf(path, width = width, height = height)
    on.exit(dev.off(), add = TRUE)
    plot_obj <- force(expr)
    if (!is.null(plot_obj)) print(plot_obj)
}

read_tsv <- function(path) {
    read.table(
        path, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
        quote = "", fill = TRUE, comment.char = "", check.names = FALSE
    )
}

tsv_update_annotation <- function(tsv, accession, annotation_path) {
    system2(
        "python3",
        args   = c("src/tsv_updater.py", tsv, accession, "Annotated=TRUE",
                   paste0("Annotation_Path=", annotation_path)),
        stdout = FALSE,
        stderr = FALSE
    )
}

prepare_peaks <- function(peak_file) {
    peaks <- readPeakFile(peak_file, as = "GRanges")
    if (length(peaks) == 0L) return(peaks)
    if (!any(grepl("^chr", GenomeInfoDb::seqlevels(peaks)))) {
        GenomeInfoDb::seqlevels(peaks) <- paste0("chr", GenomeInfoDb::seqlevels(peaks))
    }
    GenomeInfoDb::keepSeqlevels(
        peaks,
        intersect(GenomeInfoDb::seqlevels(peaks), main_chrs),
        pruning.mode = "coarse"
    )
}

plot_avg_profile <- function(peaks, label, sample_fig_dir) {
    gr <- peaks

    if (length(gr) > max_peaks_prof) {
        gr <- gr[seq_len(max_peaks_prof)]
        message(sprintf("[avgprof] %s — truncated to %d peaks", label, max_peaks_prof))
    }

    if (length(gr) == 0L) {
        message(sprintf("[avgprof] %s — no peaks, skipping", label))
        return(invisible(FALSE))
    }

    tag_matrix <- tryCatch(
        getTagMatrix(gr, windows = promoter_windows),
        error = function(e) {
            message(sprintf("[avgprof] %s — getTagMatrix failed: %s", label, conditionMessage(e)))
            NULL
        }
    )

    if (is.null(tag_matrix)) return(invisible(FALSE))

    tag_matrix <- tag_matrix[rowSums(tag_matrix) > 0, , drop = FALSE]

    if (nrow(tag_matrix) < 2L) {
        message(sprintf("[avgprof] %s — fewer than 2 non-zero rows, skipping", label))
        return(invisible(FALSE))
    }

    tryCatch(
        save_pdf(
            plotAvgProf(tag_matrix, xlim = c(-tss_upstream, tss_downstream),
                        conf = FALSE, xlab = "Distance to TSS (bp)", ylab = "Read count frequency"),
            file.path(sample_fig_dir, paste0(label, "_avgprof.pdf")),
            width = 8, height = 5
        ),
        error = function(e) message(sprintf("[avgprof] %s — plotAvgProf failed: %s", label, conditionMessage(e)))
    )

    rm(tag_matrix); gc()
    invisible(TRUE)
}

process_sample <- function(label, rows, tsv_path, anno_root, fig_root) {
    peak_file       <- rows$Peak_File[1]
    accessions      <- unique(rows$Accession)
    sample_anno_dir <- file.path(anno_root, label)
    sample_fig_dir  <- file.path(fig_root, label)
    annotation_path <- file.path(sample_anno_dir, paste0(label, "_annotation.tsv"))
    rdata_path      <- file.path(sample_anno_dir, paste0(label, "_anno.RData"))

    if (file.exists(annotation_path) && file.exists(rdata_path)) {
        message(sprintf("[skip] %s — already annotated", label))
        for (acc in accessions) tsv_update_annotation(tsv_path, acc, annotation_path)
        return(label)
    }

    if (!file.exists(peak_file)) {
        message(sprintf("[warn] Peak file not found for %s — skipping: %s", label, peak_file))
        return(NULL)
    }

    message(sprintf("[annotate] %s", label))
    dir.create(sample_anno_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(sample_fig_dir,  recursive = TRUE, showWarnings = FALSE)

    peaks <- prepare_peaks(peak_file)
    if (length(peaks) == 0L) {
        message(sprintf("[warn] %s — no peaks on main chromosomes", label))
        return(NULL)
    }

    peak_anno <- suppressWarnings(
        annotatePeak(peaks, tssRegion = c(-tss_upstream, tss_downstream),
                     TxDb = txdb, annoDb = "org.Hs.eg.db", verbose = FALSE)
    )

    anno_df           <- as.data.frame(peak_anno)
    anno_df$sample_id <- label
    anno_df$peak_file <- peak_file

    write.table(anno_df, annotation_path, sep = "\t", quote = FALSE, row.names = FALSE)

    for (item in list(
        list(fn = plotAnnoPie,   suffix = "_annopie.pdf"),
        list(fn = plotAnnoBar,   suffix = "_annobar.pdf"),
        list(fn = plotDistToTSS, suffix = "_dist_to_tss.pdf")
    )) {
        tryCatch(
            save_pdf(item$fn(peak_anno), file.path(sample_fig_dir, paste0(label, item$suffix))),
            error = function(e) message(sprintf("[plot] %s — %s failed: %s", label, item$suffix, conditionMessage(e)))
        )
    }

    plot_avg_profile(peaks, label, sample_fig_dir)
    save(peak_anno, file = rdata_path)

    for (acc in accessions) tsv_update_annotation(tsv_path, acc, annotation_path)
    message(sprintf("[done] %s — %d peaks", label, length(peaks)))

    rm(peaks, peak_anno, anno_df); gc()
    label
}

all_processed_labels <- list()

for (tsv_path in tsv_paths) {
    if (!file.exists(tsv_path)) {
        message(sprintf("[warn] TSV not found: %s — skipping", tsv_path))
        next
    }

    tsv_base  <- sub("\\.tsv$", "", basename(tsv_path))
    anno_root <- file.path("pipeline_outputs", "chipseq", tsv_base, "annotation")
    fig_root  <- file.path("results", "figures", "chipseq", tsv_base, "annotation")
    tab_root  <- file.path("results", "tables",  "chipseq", tsv_base, "annotation")

    dir.create(anno_root, recursive = TRUE, showWarnings = FALSE)
    dir.create(fig_root,  recursive = TRUE, showWarnings = FALSE)
    dir.create(tab_root,  recursive = TRUE, showWarnings = FALSE)

    meta <- read_tsv(tsv_path)

    required_cols <- c("Accession", "Condition", "Sample_Type", "Replicate", "Peak_File")
    missing_cols  <- setdiff(required_cols, colnames(meta))
    if (length(missing_cols) > 0L) {
        message(sprintf("[warn] %s — missing columns: %s — skipping",
                        tsv_base, paste(missing_cols, collapse = ", ")))
        next
    }

    samples <- meta[tolower(meta$Sample_Type) == "ip" &
                    !is.na(meta$Peak_File) & nzchar(meta$Peak_File), ]

    if (nrow(samples) == 0L) {
        message(sprintf("[warn] %s — no IP rows with Peak_File, skipping", tsv_base))
        next
    }

    samples$sample_id <- mapply(sample_id, samples$Condition, samples$Sample_Type,
                                samples$Replicate, USE.NAMES = FALSE)

    message(sprintf("[setup] %s | %d IP row(s), %d sample(s)",
                    tsv_base, nrow(samples), length(unique(samples$sample_id))))

    processed_labels <- c()

    for (label in unique(samples$sample_id)) {
        rows <- samples[samples$sample_id == label, ]
        if (length(unique(rows$Peak_File)) > 1L) {
            message(sprintf("[warn] %s has multiple Peak_File values; using first", label))
            rows <- rows[rows$Peak_File == rows$Peak_File[1], ]
        }
        processed <- process_sample(label, rows, tsv_path, anno_root, fig_root)
        if (!is.null(processed)) processed_labels <- c(processed_labels, processed)
        gc()
    }

    if (length(processed_labels) == 0L) {
        message(sprintf("[warn] %s — no samples successfully annotated", tsv_base))
        next
    }

    if (length(processed_labels) > 1L) {
        anno_list <- list()
        for (label in processed_labels) {
            rdata_path <- file.path(anno_root, label, paste0(label, "_anno.RData"))
            if (file.exists(rdata_path)) {
                local({ load(rdata_path); anno_list[[label]] <<- peak_anno })
            }
            gc()
        }
        tryCatch(
            save_pdf(plotAnnoBar(anno_list), file.path(fig_root, "annobar_comparative.pdf"),
                     width = 9, height = 5),
            error = function(e) message(sprintf("[comparative] plotAnnoBar failed: %s", conditionMessage(e)))
        )
        rm(anno_list); gc()
    }

    annotation_tables <- Filter(Negate(is.null), lapply(processed_labels, function(label) {
        f <- file.path(anno_root, label, paste0(label, "_annotation.tsv"))
        if (file.exists(f)) read_tsv(f) else NULL
    }))

    if (length(annotation_tables) > 0L) {
        all_anno_df <- do.call(rbind, annotation_tables)
        rm(annotation_tables); gc()

        for (out_path in c(file.path(anno_root, "all_annotations.tsv"),
                           file.path(tab_root,  "all_annotations.tsv"))) {
            write.table(all_anno_df, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
        }
        message(sprintf("[export] Combined annotations → %s", file.path(anno_root, "all_annotations.tsv")))

        if ("annotation" %in% colnames(all_anno_df) && "geneId" %in% colnames(all_anno_df)) {
            promoter_df    <- all_anno_df[grepl("Promoter", all_anno_df$annotation), , drop = FALSE]
            promoter_genes <- unique(as.character(promoter_df$geneId))
            promoter_genes <- promoter_genes[!is.na(promoter_genes) & nzchar(promoter_genes)]
        } else {
            promoter_df    <- data.frame()
            promoter_genes <- character(0)
        }

        for (out_path in c(file.path(anno_root, "promoter_genes.txt"),
                           file.path(tab_root,  "promoter_genes.txt"))) {
            writeLines(promoter_genes, out_path)
        }
        message(sprintf("[export] %d unique promoter genes → %s",
                        length(promoter_genes), file.path(anno_root, "promoter_genes.txt")))

        if (nrow(promoter_df) > 0L) {
            peak_counts <- sort(table(promoter_df$geneId), decreasing = TRUE)
            multi_peak  <- peak_counts[peak_counts > 1]
            if (length(multi_peak) > 0L) {
                multi_peak_df <- data.frame(geneId = names(multi_peak),
                                            promoter_peak_count = as.integer(multi_peak),
                                            row.names = NULL)
                for (out_path in c(file.path(anno_root, "multi_peak_promoter_genes.tsv"),
                                   file.path(tab_root,  "multi_peak_promoter_genes.tsv"))) {
                    write.table(multi_peak_df, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
                }
                message(sprintf("[multi-peak] %d gene(s) with >1 promoter peak", nrow(multi_peak_df)))
            }
        }

        rm(all_anno_df, promoter_df); gc()
    }

    all_processed_labels[[tsv_base]] <- processed_labels
}

total <- sum(lengths(all_processed_labels))
if (total == 0L) stop("No samples were successfully annotated.", call. = FALSE)

message(sprintf("[done] Peak annotation pipeline complete. %d sample(s) processed.", total))