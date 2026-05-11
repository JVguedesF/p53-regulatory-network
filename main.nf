nextflow.enable.dsl = 2

include { DOWNLOAD_METADATA; DOWNLOAD_GENOME; DOWNLOAD_CDNA } from './modules/setup'
include { QC_TRIMMING } from './modules/preprocessing'
include { CHIP_ALIGNMENT; DEDUPLICATION; PEAK_CALLING; PEAK_ANNOTATION } from './modules/chipseq'
include { RNA_QUANTIFICATION; DESEQ2 } from './modules/rnaseq'
include { INTEGRATION; MULTILINHAGEM; SURVIVAL_TCGA; IMMUNE_TARGETS; PUBLICATION_FIGURES } from './modules/downstream'

def classifyTsv(Path p) {
  def lines = p.readLines()

  if (lines.size() < 2) {
    return [false, false]
  }

  def header = lines[0].split('\t', -1) as List
  def idx = header.indexOf("Sample_Type")

  if (idx < 0) {
    return [false, false]
  }

  def types = lines.drop(1).collect { line ->
    def cols = line.split('\t', -1)
    cols.size() > idx ? cols[idx].toLowerCase() : ""
  }

  def hasRna  = types.any { sampleType -> sampleType == "rna" }
  def hasChip = types.any { sampleType -> sampleType == "ip" || sampleType == "input" }

  return [hasRna, hasChip]
}

workflow {
  repo_scripts = file("${projectDir}/scripts", checkIfExists: true)
  repo_src     = file("${projectDir}/src", checkIfExists: true)
  repo_config  = file("${projectDir}/config", checkIfExists: true)
  project_root = projectDir.toString()

  if (params.run_download) {
    downloaded = DOWNLOAD_METADATA(repo_scripts, repo_src, repo_config, project_root)

    tsv_ch = downloaded.tsv_list
      .splitText()
      .map { line ->
        def clean = line.trim()
        def tsv_file = file(clean)
        def classes = classifyTsv(tsv_file)
        tuple(tsv_file.toRealPath().toString(), tsv_file.baseName, classes[0], classes[1])
      }
  } else {
    tsv_ch = channel
      .fromPath("${params.samples_tsv}/*.tsv", checkIfExists: true)
      .map { tsv_file ->
        def classes = classifyTsv(tsv_file)
        tuple(tsv_file.toRealPath().toString(), tsv_file.baseName, classes[0], classes[1])
      }
  }

  genome = DOWNLOAD_GENOME(repo_scripts, repo_src, repo_config, project_root)
  cdna   = DOWNLOAD_CDNA(repo_scripts, repo_src, repo_config, project_root)

  trimmed = QC_TRIMMING(tsv_ch, repo_scripts, repo_src, repo_config, project_root)

  chip_trimmed_ch = trimmed.done
    .filter { _tsv, _dataset, _has_rna, has_chip, _done -> has_chip }
    .map { tsv, dataset, _has_rna, _has_chip, _done -> tuple(tsv, dataset) }

  rna_trimmed_ch = trimmed.done
    .filter { _tsv, _dataset, has_rna, _has_chip, _done -> has_rna }
    .map { tsv, dataset, _has_rna, _has_chip, _done -> tuple(tsv, dataset) }

  chip_ready_ch = chip_trimmed_ch.combine(genome.done)
  rna_ready_ch  = rna_trimmed_ch.combine(cdna.done)

  aligned   = CHIP_ALIGNMENT(chip_ready_ch, repo_scripts, repo_src, repo_config, project_root)
  deduped   = DEDUPLICATION(aligned.done, repo_scripts, repo_src, repo_config, project_root)
  peaks     = PEAK_CALLING(deduped.done, repo_scripts, repo_src, repo_config, project_root)
  annotated = PEAK_ANNOTATION(peaks.done, repo_scripts, repo_src, repo_config, project_root)

  quantified = RNA_QUANTIFICATION(rna_ready_ch, repo_scripts, repo_src, repo_config, project_root)
  deseq2     = DESEQ2(quantified.done, repo_scripts, repo_src, repo_config, project_root)

  if (params.run_downstream) {
    peak_annotation_done_files = annotated.done
      .map { _tsv, _dataset, done -> done }
      .collect()

    deseq2_done_files = deseq2.done
      .map { _tsv, _dataset, done -> done }
      .collect()

    pairs_csv = channel.fromPath(params.integration_pairs, checkIfExists: true)

    integration = INTEGRATION(
      pairs_csv,
      peak_annotation_done_files,
      deseq2_done_files,
      repo_scripts,
      repo_src,
      repo_config,
      project_root
    )

    multilinhagem = MULTILINHAGEM(
      integration.done,
      repo_scripts,
      repo_src,
      repo_config,
      project_root
    )

    survival = SURVIVAL_TCGA(
      multilinhagem.done,
      repo_scripts,
      repo_src,
      repo_config,
      project_root
    )

    immune = IMMUNE_TARGETS(
      multilinhagem.done,
      survival.done,
      repo_scripts,
      repo_src,
      repo_config,
      project_root
    )

    PUBLICATION_FIGURES(
      integration.done,
      multilinhagem.done,
      survival.done,
      immune.done,
      repo_scripts,
      repo_src,
      repo_config,
      project_root
    )
  }
}