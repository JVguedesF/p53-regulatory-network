process DOWNLOAD_METADATA {
  label "low_mem"

  input:
  path scripts_dir
  path src_dir
  path config_dir
  val project_root

  output:
  path "tsv_files.txt", emit: tsv_list

  script:
  """
  set -euo pipefail

  WORKDIR="\$PWD"
  : "${scripts_dir}" "${src_dir}" "${config_dir}"
  cd "${project_root}"

  mkdir -p "${params.samples_tsv}" "${params.logs}"

  export PIPELINE_OUTPUTS="${params.pipeline_outputs}"
  export RESULTS="${params.results}"
  export LOGS="${params.logs}"
  export OUT_DIR="${params.pipeline_outputs}"
  export RESULTS_DIR="${params.results}"
  export LOG_DIR="${params.logs}"

  python3 scripts/00_download_data.py

  find "${params.samples_tsv}" -maxdepth 1 -type f -name "*.tsv" | sort > "\$WORKDIR/tsv_files.txt"
  """
}

process DOWNLOAD_GENOME {
  label "low_mem"

  input:
  path scripts_dir
  path src_dir
  path config_dir
  val project_root

  output:
  path "genome.done", emit: done

  script:
  """
  set -euo pipefail

  WORKDIR="\$PWD"
  : "${scripts_dir}" "${src_dir}" "${config_dir}"
  cd "${project_root}"

  mkdir -p "\$(dirname "${params.genome_fasta}")"

  export PIPELINE_OUTPUTS="${params.pipeline_outputs}"
  export RESULTS="${params.results}"
  export LOGS="${params.logs}"
  export OUT_DIR="${params.pipeline_outputs}"
  export RESULTS_DIR="${params.results}"
  export LOG_DIR="${params.logs}"

  python3 scripts/01_download_genome.py

  test -s "${params.genome_fasta}"

  touch "\$WORKDIR/genome.done"
  """
}

process DOWNLOAD_CDNA {
  label "low_mem"

  input:
  path scripts_dir
  path src_dir
  path config_dir
  val project_root

  output:
  path "cdna.done", emit: done

  script:
  """
  set -euo pipefail

  WORKDIR="\$PWD"
  : "${scripts_dir}" "${src_dir}" "${config_dir}"
  cd "${project_root}"

  mkdir -p "\$(dirname "${params.cdna_fasta}")"

  export PIPELINE_OUTPUTS="${params.pipeline_outputs}"
  export RESULTS="${params.results}"
  export LOGS="${params.logs}"
  export OUT_DIR="${params.pipeline_outputs}"
  export RESULTS_DIR="${params.results}"
  export LOG_DIR="${params.logs}"

  if [[ ! -s "${params.cdna_fasta}" ]]; then
    python3 - <<'PY'
import os
import urllib.request

url = "${params.cdna_url}"
out = "${params.cdna_fasta}"

os.makedirs(os.path.dirname(out), exist_ok=True)
urllib.request.urlretrieve(url, out)
PY
  fi

  test -s "${params.cdna_fasta}"

  touch "\$WORKDIR/cdna.done"
  """
}
