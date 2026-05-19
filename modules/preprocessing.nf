process QC_TRIMMING {
  tag "$dataset"
  label "low_mem"

  input:
  tuple val(tsv), val(dataset), val(has_rna), val(has_chip)
  path scripts_dir
  path src_dir
  path config_dir
  val project_root

  output:
  tuple val(tsv), val(dataset), val(has_rna), val(has_chip), path("${dataset}.qc.done"), emit: done

  script:
  """
  set -euo pipefail

  WORKDIR="\$PWD"
  : "${scripts_dir}" "${src_dir}" "${config_dir}"
  cd "${project_root}"

  mkdir -p "${params.pipeline_outputs}" "${params.results}" "${params.logs}"

  export PIPELINE_OUTPUTS="${params.pipeline_outputs}"
  export RESULTS="${params.results}"
  export LOGS="${params.logs}"
  export OUT_DIR="${params.pipeline_outputs}"
  export RESULTS_DIR="${params.results}"
  export LOG_DIR="${params.logs}"

  bash scripts/02_qc_trimming.sh "${tsv}"

  touch "\$WORKDIR/${dataset}.qc.done"
  """
}
