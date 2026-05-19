process CHIP_ALIGNMENT {
  tag "$dataset"
  label "high_mem"

  input:
  tuple val(tsv), val(dataset), path(genome_done)
  path scripts_dir
  path src_dir
  path config_dir
  val project_root

  output:
  tuple val(tsv), val(dataset), path("${dataset}.alignment.done"), emit: done

  script:
  """
  set -euo pipefail

  WORKDIR="\$PWD"
  : "${scripts_dir}" "${src_dir}" "${config_dir}"
  cd "${project_root}"

  export PIPELINE_OUTPUTS="${params.pipeline_outputs}"
  export RESULTS="${params.results}"
  export LOGS="${params.logs}"
  export OUT_DIR="${params.pipeline_outputs}"
  export RESULTS_DIR="${params.results}"
  export LOG_DIR="${params.logs}"

  bash scripts/03a_alignment_chipseq.sh "${tsv}"

  touch "\$WORKDIR/${dataset}.alignment.done"
  """
}

process DEDUPLICATION {
  tag "$dataset"
  label "high_mem"

  input:
  tuple val(tsv), val(dataset), path(alignment_done)
  path scripts_dir
  path src_dir
  path config_dir
  val project_root

  output:
  tuple val(tsv), val(dataset), path("${dataset}.dedup.done"), emit: done

  script:
  """
  set -euo pipefail

  WORKDIR="\$PWD"
  : "${scripts_dir}" "${src_dir}" "${config_dir}"
  cd "${project_root}"

  export PIPELINE_OUTPUTS="${params.pipeline_outputs}"
  export RESULTS="${params.results}"
  export LOGS="${params.logs}"
  export OUT_DIR="${params.pipeline_outputs}"
  export RESULTS_DIR="${params.results}"
  export LOG_DIR="${params.logs}"

  bash scripts/04_deduplication.sh "${tsv}"

  touch "\$WORKDIR/${dataset}.dedup.done"
  """
}

process PEAK_CALLING {
  tag "$dataset"
  label "medium_mem"

  input:
  tuple val(tsv), val(dataset), path(dedup_done)
  path scripts_dir
  path src_dir
  path config_dir
  val project_root

  output:
  tuple val(tsv), val(dataset), path("${dataset}.peaks.done"), emit: done

  script:
  """
  set -euo pipefail

  WORKDIR="\$PWD"
  : "${scripts_dir}" "${src_dir}" "${config_dir}"
  cd "${project_root}"

  export PIPELINE_OUTPUTS="${params.pipeline_outputs}"
  export RESULTS="${params.results}"
  export LOGS="${params.logs}"
  export OUT_DIR="${params.pipeline_outputs}"
  export RESULTS_DIR="${params.results}"
  export LOG_DIR="${params.logs}"

  bash scripts/05_peak_calling.sh "${tsv}"

  touch "\$WORKDIR/${dataset}.peaks.done"
  """
}

process PEAK_ANNOTATION {
  tag "$dataset"
  label "medium_mem"

  input:
  tuple val(tsv), val(dataset), path(peaks_done)
  path scripts_dir
  path src_dir
  path config_dir
  val project_root

  output:
  tuple val(tsv), val(dataset), path("${dataset}.peak_annotation.done"), emit: done

  script:
  """
  set -euo pipefail

  WORKDIR="\$PWD"
  : "${scripts_dir}" "${src_dir}" "${config_dir}"
  cd "${project_root}"

  export PIPELINE_OUTPUTS="${params.pipeline_outputs}"
  export RESULTS="${params.results}"
  export LOGS="${params.logs}"
  export OUT_DIR="${params.pipeline_outputs}"
  export RESULTS_DIR="${params.results}"
  export LOG_DIR="${params.logs}"

  Rscript scripts/07_peak_annotation.R "${tsv}"

  touch "\$WORKDIR/${dataset}.peak_annotation.done"
  """
}
