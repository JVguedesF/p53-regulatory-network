nextflow.enable.dsl=2

process DOWNLOAD {
  label "low_mem"
  publishDir "${params.outdir}", mode: "copy"

  input:
  path tsv

  output:
  path tsv, emit: tsv

  script:
  """
  python scripts/00_download_data.py
  """
}

process DOWNLOAD_GENOME {
  label "low_mem"

  output:
  path "data/genome/fasta/hg38.fa", emit: fasta

  script:
  """
  python scripts/01_download_genome.py
  """
}

process QC_TRIMMING {
  tag "$tsv.simpleName"
  label "low_mem"

  input:
  path tsv

  output:
  path tsv, emit: tsv

  script:
  """
  bash scripts/02_qc_trimming.sh ${tsv}
  """
}

process ALIGNMENT {
  tag "$tsv.simpleName"
  label "high_mem"

  input:
  path tsv
  path fasta

  output:
  path tsv, emit: tsv

  script:
  """
  bash scripts/03_alignment.sh ${tsv}
  """
}

process DEDUPLICATION {
  tag "$tsv.simpleName"
  label "high_mem"

  input:
  path tsv

  output:
  path tsv, emit: tsv

  script:
  """
  bash scripts/04_deduplication.sh ${tsv}
  """
}

process PEAK_CALLING {
  tag "$tsv.simpleName"
  label "low_mem"

  input:
  path tsv

  output:
  path tsv, emit: tsv

  script:
  """
  bash scripts/05_peak_calling.sh ${tsv}
  """
}

process PEAK_ANNOTATION {
  tag "$tsv.simpleName"
  label "low_mem"
  publishDir "${params.outdir}/annotation", mode: "copy"

  input:
  path tsv

  output:
  path "results/annotation/**", emit: annotation
  path "results/figures/**",    emit: figures
  path "results/peaks/**",      emit: gene_lists

  script:
  """
  Rscript scripts/06_peak_annotation.R ${tsv}
  """
}

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------

workflow {
  tsv_ch = channel
    .fromPath("${params.samples_tsv}/*.tsv", checkIfExists: true)

  genome     = DOWNLOAD_GENOME()
  downloaded = DOWNLOAD(tsv_ch)

  trimmed = QC_TRIMMING(downloaded.tsv)
  aligned = ALIGNMENT(trimmed.tsv, genome.fasta)
  deduped = DEDUPLICATION(aligned.tsv)
  peaks   = PEAK_CALLING(deduped.tsv)

  PEAK_ANNOTATION(peaks.tsv)
}