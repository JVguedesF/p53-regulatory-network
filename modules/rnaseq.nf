process RNA_QUANTIFICATION {
  tag "$dataset"
  label "high_mem"

  input:
  tuple val(tsv), val(dataset), path(cdna_done)
  path scripts_dir
  path src_dir
  path config_dir
  val project_root

  output:
  tuple val(tsv), val(dataset), path("${dataset}.rna_quant.done"), emit: done

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

  bash scripts/03b_quantification_rnaseq.sh "${tsv}"

  touch "\$WORKDIR/${dataset}.rna_quant.done"
  """
}

process DESEQ2 {
  tag "$dataset"
  label "high_mem"

  input:
  tuple val(tsv), val(dataset), path(rna_quant_done)
  path scripts_dir
  path src_dir
  path config_dir
  val project_root

  output:
  tuple val(tsv), val(dataset), path("${dataset}.deseq2.done"), emit: done

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

  if [[ -n "${params.deseq_treated}" && -n "${params.deseq_reference}" ]]; then
    Rscript scripts/08_differential_expression.R \
      "${tsv}" \
      "${params.deseq_treated}" \
      "${params.deseq_reference}" \
      "${params.lfc_threshold}" \
      "${params.padj_threshold}"
  else
    Rscript scripts/08_differential_expression.R "${tsv}"
  fi

  touch "\$WORKDIR/${dataset}.deseq2.done"
  """
}
