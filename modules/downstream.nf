process INTEGRATION {
  tag "integration"
  label "medium_mem"

  input:
  path pairs_csv
  path peak_annotation_done_files
  path deseq2_done_files
  path scripts_dir
  path src_dir
  path config_dir
  val project_root

  output:
  path "integration.done", emit: done

  script:
  """
  set -euo pipefail

  PAIRS_CSV="\$(realpath "${pairs_csv}")"

  WORKDIR="\$PWD"
  : "${scripts_dir}" "${src_dir}" "${config_dir}"
  cd "${project_root}"

  export PIPELINE_OUTPUTS="${params.pipeline_outputs}"
  export RESULTS="${params.results}"
  export LOGS="${params.logs}"
  export OUT_DIR="${params.pipeline_outputs}"
  export RESULTS_DIR="${params.results}"
  export LOG_DIR="${params.logs}"
  export PAIRS_CSV
  export PROJECT_DIR="${project_root}"
  export SAMPLES_TSV="${params.samples_tsv}"
  export LFC_THRESHOLD="${params.lfc_threshold}"
  export PADJ_THRESHOLD="${params.padj_threshold}"
  export PROMOTER_WINDOW="${params.promoter_window}"

  python3 - <<'PY'
import csv
import os
import subprocess
from pathlib import Path

project_dir = Path(os.environ["PROJECT_DIR"])
samples_tsv = Path(os.environ["SAMPLES_TSV"])
pairs_csv   = Path(os.environ["PAIRS_CSV"])

lfc      = os.environ["LFC_THRESHOLD"]
padj     = os.environ["PADJ_THRESHOLD"]
promoter = os.environ["PROMOTER_WINDOW"]

def get(row, *names):
    for name in names:
        value = (row.get(name) or "").strip()
        if value:
            return value
    return ""

def resolve_tsv(row, explicit_field, dataset_field):
    explicit = get(row, explicit_field)
    if explicit:
        p = Path(explicit)
        if not p.is_absolute():
            p = project_dir / p
        return p.resolve()

    dataset = get(row, dataset_field)
    if not dataset:
        raise ValueError(f"Missing {explicit_field} or {dataset_field} in integration_pairs.csv")

    return (samples_tsv / f"{dataset}.tsv").resolve()

with pairs_csv.open(newline="") as fh:
    reader = csv.DictReader(fh)

    for row in reader:
        chip_tsv = resolve_tsv(row, "chip_tsv", "chip_dataset")
        rna_tsv  = resolve_tsv(row, "rnaseq_tsv", "rna_dataset")

        chip_dataset = get(row, "chip_dataset") or chip_tsv.stem
        rna_dataset  = get(row, "rna_dataset")  or rna_tsv.stem

        cmd = [
            "Rscript",
            "scripts/09_integration.R",
            str(chip_tsv),
            str(rna_tsv),
            chip_dataset,
            rna_dataset,
            lfc,
            padj,
            promoter,
        ]

        print("[integration]", " ".join(cmd), flush=True)
        subprocess.run(cmd, check=True)
PY

  touch "\$WORKDIR/integration.done"
  """
}

process MULTILINHAGEM {
  tag "multilinhagem"
  label "medium_mem"

  input:
  path integration_done
  path scripts_dir
  path src_dir
  path config_dir
  val project_root

  output:
  path "multilinhagem.done", emit: done

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

  Rscript scripts/10_multilinhagem.R \
    "${params.pipeline_outputs}/integration" \
    "${params.lfc_threshold}" \
    "${params.padj_threshold}"

  touch "\$WORKDIR/multilinhagem.done"
  """
}

process SURVIVAL_TCGA {
  tag "$params.cancer_type"
  label "high_mem"

  input:
  path multilinhagem_done
  path scripts_dir
  path src_dir
  path config_dir
  val project_root

  output:
  path "survival.done", emit: done

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

  Rscript scripts/11_survival_tcga.R "${params.cancer_type}"

  touch "\$WORKDIR/survival.done"
  """
}

process IMMUNE_TARGETS {
  tag "immune_targets"
  label "low_mem"

  input:
  path multilinhagem_done
  path survival_done
  path scripts_dir
  path src_dir
  path config_dir
  val project_root

  output:
  path "immune_targets.done", emit: done

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

  python3 scripts/12_immune_targets.py \
    --conserved "${params.pipeline_outputs}/multilinhagem/conserved_targets_all_runs.tsv" \
    --survival "${params.results}/tables/survival/${params.cancer_type}/cox_results.tsv" \
    --cancer-type "${params.cancer_type}" \
    --fdr-survival "${params.padj_threshold}"

  touch "\$WORKDIR/immune_targets.done"
  """
}

process PUBLICATION_FIGURES {
  tag "publication"
  label "medium_mem"

  input:
  path integration_done
  path multilinhagem_done
  path survival_done
  path immune_targets_done
  path scripts_dir
  path src_dir
  path config_dir
  val project_root

  output:
  path "publication_figures.done", emit: done

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

  Rscript scripts/13_figures.R "${params.cancer_type}"

  touch "\$WORKDIR/publication_figures.done"
  """
}
