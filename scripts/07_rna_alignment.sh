#!/bin/bash
set -euo pipefail
shopt -s nullglob

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$(dirname "$SCRIPT_DIR")"

(( BASH_VERSINFO[0] >= 4 )) || { echo "Bash 4+ required (found $BASH_VERSION)" >&2; exit 1; }

TSV_DIR="data/tsv"
OUT_DIR="pipeline_outputs"
LOG_DIR="logs"
TARGET="hg38"
CDNA_FASTA="data/genome/fasta/Homo_sapiens.GRCh38.cdna.all.fa.gz"
SALMON_INDEX="data/genome/index/salmon/salmon_index_${TARGET}"
THREADS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
MAX_JOBS=1

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; NC='\033[0m'

log()     { echo "[$(date '+%H:%M:%S')] $*"; [[ -n "${LOG:-}" ]] && echo "[$(date '+%H:%M:%S')] $*" >> "$LOG"; return 0; }
die()     { echo -e "${RED}FATAL ERROR:${NC} $*" >&2; [[ -n "${LOG:-}" ]] && echo "FATAL ERROR: $*" >> "$LOG"; exit 1; }
warn()    { echo -e "${YELLOW}WARNING:${NC} $*" >&2; [[ -n "${LOG:-}" ]] && echo "WARNING: $*" >> "$LOG"; return 0; }
success() { echo -e "${GREEN}OK${NC} $*"; [[ -n "${LOG:-}" ]] && echo "OK $*" >> "$LOG"; return 0; }

tsv_update() {
    local tsv="$1" acc="$2"; shift 2
    flock -x "${tsv}.lock" python src/tsv_updater.py "$tsv" "$acc" "$@"
}

write_log_header() {
    local mode="$1" id="$2"
    { flock 200
      printf '\n%s\n' '====================================================================='
      echo "[$(date '+%H:%M:%S')] START ${mode}: ${id}"
      printf '%s\n' '---------------------------------------------------------------------'
    } >> "$LOG" 200>> "$LOG"
}

quant_se() {
    local id="$1" acc="$2" trimmed="$3" tsv="$4"
    local safe_id="${id// /_}"
    local quant_dir="${QUANT_OUT}/${safe_id}_${acc}"

    if [[ -f "${quant_dir}/quant.sf" ]]; then
        log "Skipping SE: $id — already quantified."
        tsv_update "$tsv" "$acc" "BAM_Path=${quant_dir}"
        return 0
    fi

    write_log_header "QUANT_SE" "$id"

    local tmp; tmp=$(mktemp)
    salmon quant \
        --index          "$SALMON_INDEX" \
        --libType        A \
        --unmatedReads   "$trimmed" \
        --threads        "$THREADS" \
        --validateMappings \
        --gcBias \
        --output         "$quant_dir" > "$tmp" 2>&1 \
    || { echo "[$(date '+%H:%M:%S')] ERROR $id" >> "$LOG"; cat "$tmp" >> "$LOG"; rm -f "$tmp"; exit 1; }

    { flock 200
      cat "$tmp"
      echo "[$(date '+%H:%M:%S')] END QUANT_SE: $id"
    } >> "$LOG" 200>> "$LOG"
    rm -f "$tmp"

    [[ -f "${quant_dir}/quant.sf" ]] || die "Salmon did not produce quant.sf for $id"

    tsv_update "$tsv" "$acc" "BAM_Path=${quant_dir}"
    success "$id — quantification done → ${quant_dir}"
}

quant_pe() {
    local id="$1" r1="$2" r2="$3" acc1="$4" acc2="$5" tsv="$6"
    local safe_id="${id// /_}"
    local acc_base="${acc1%_1}"
    local quant_dir="${QUANT_OUT}/${safe_id}_${acc_base}"

    if [[ -f "${quant_dir}/quant.sf" ]]; then
        log "Skipping PE: $id — already quantified."
        tsv_update "$tsv" "$acc1" "BAM_Path=${quant_dir}"
        tsv_update "$tsv" "$acc2" "BAM_Path=${quant_dir}"
        return 0
    fi

    write_log_header "QUANT_PE" "$id"

    local tmp; tmp=$(mktemp)
    salmon quant \
        --index          "$SALMON_INDEX" \
        --libType        A \
        --mates1         "$r1" \
        --mates2         "$r2" \
        --threads        "$THREADS" \
        --validateMappings \
        --gcBias \
        --output         "$quant_dir" > "$tmp" 2>&1 \
    || { echo "[$(date '+%H:%M:%S')] ERROR $id" >> "$LOG"; cat "$tmp" >> "$LOG"; rm -f "$tmp"; exit 1; }

    { flock 200
      cat "$tmp"
      echo "[$(date '+%H:%M:%S')] END QUANT_PE: $id"
    } >> "$LOG" 200>> "$LOG"
    rm -f "$tmp"

    [[ -f "${quant_dir}/quant.sf" ]] || die "Salmon did not produce quant.sf for $id"

    tsv_update "$tsv" "$acc1" "BAM_Path=${quant_dir}"
    tsv_update "$tsv" "$acc2" "BAM_Path=${quant_dir}"
    success "$id — quantification done → ${quant_dir}"
}

for cmd in salmon python; do
    command -v "$cmd" &>/dev/null || die "'$cmd' not found in PATH."
done
python -c "import src.tsv_updater" 2>/dev/null \
    || die "src/tsv_updater.py not importable — run from project root."

mkdir -p "$(dirname "$SALMON_INDEX")" "$LOG_DIR"

if [[ ! -d "$SALMON_INDEX" ]]; then
    [[ -f "$CDNA_FASTA" ]] || die \
        "cDNA FASTA not found: ${CDNA_FASTA}
        Download with:
          wget -P data/genome/fasta/ \\
            'https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz'"

    log "Building Salmon index for $TARGET..."
    salmon index \
        --transcripts "$CDNA_FASTA" \
        --index       "$SALMON_INDEX" \
        --threads     "$THREADS" \
        2>&1 | tee -a "$LOG_DIR/salmon_index.log" \
        || die "Failed to build Salmon index."
    success "Salmon index ready → $SALMON_INDEX"
else
    log "Salmon index found — skipping build."
fi

if [[ $# -gt 0 ]]; then
    TSV_FILES=("$@")
else
    mapfile -t TSV_FILES < <(find "$TSV_DIR" -maxdepth 1 -name "*.tsv" | sort)
fi
[[ ${#TSV_FILES[@]} -eq 0 ]] && die "No .tsv files found in $TSV_DIR"

for TSV_FILE in "${TSV_FILES[@]}"; do
    [[ -f "$TSV_FILE" ]] || { warn "$TSV_FILE not found — skipping."; continue; }
    TSV_BASE="${TSV_FILE##*/}"; TSV_BASE="${TSV_BASE%.tsv}"

    has_rna=false
    declare -a se_order=()
    declare -A se_acc=() se_trimmed=()
    declare -A pe_r1=() pe_r2=() pe_acc1=() pe_acc2=()

    line_num=1
    while IFS=$'\x01' read -r cond acc sample_type pair_id paired_end rep \
                                 _input_acc _fastq _dl _sz _md5 \
                                 _qc trimmed_path _rest; do
        (( line_num++ ))
        [[ -z "$acc" || -z "$sample_type" ]] && continue
        [[ "${sample_type,,}" != "rna" ]] && continue

        has_rna=true

        [[ -f "$trimmed_path" ]] \
            || die "Line $line_num: Trimmed_Path not found → '$trimmed_path' (run script 02 first?)"

        group="${cond}_RNA_Rep${rep}"
        if [[ "${paired_end,,}" == "true" ]]; then
            case "$pair_id" in
                1) pe_r1["$group"]="$trimmed_path"; pe_acc1["$group"]="$acc" ;;
                2) pe_r2["$group"]="$trimmed_path"; pe_acc2["$group"]="$acc" ;;
                *) die "Line $line_num: unexpected Pair_ID='$pair_id'." ;;
            esac
        else
            se_order+=("$group")
            se_acc["$group"]="$acc"
            se_trimmed["$group"]="$trimmed_path"
        fi
    done < <(tail -n +2 "$TSV_FILE" | tr $'\t' $'\x01')

    if ! "$has_rna"; then
        log "No RNA rows in $TSV_FILE — skipping."
        continue
    fi

    QUANT_OUT="${OUT_DIR}/rnaseq/${TSV_BASE}/07_quantification"
    LOG="${LOG_DIR}/rnaseq/${TSV_BASE}_quantification.log"
    mkdir -p "$QUANT_OUT" "$(dirname "$LOG")"

    log "Processing: $TSV_FILE"

    declare -a _pids=() _ids=()
    for group in "${se_order[@]}"; do
        quant_se "$group" "${se_acc[$group]}" "${se_trimmed[$group]}" "$TSV_FILE" &
        _pids+=($!) _ids+=("$group")
        if (( ${#_pids[@]} >= MAX_JOBS )); then
            wait "${_pids[0]}" || die "SE quant failed: ${_ids[0]}"
            _pids=("${_pids[@]:1}") _ids=("${_ids[@]:1}")
        fi
    done
    for i in "${!_pids[@]}"; do wait "${_pids[$i]}" || die "SE quant failed: ${_ids[$i]}"; done

    for group in "${!pe_r1[@]}"; do
        [[ -z "${pe_r2[$group]:-}" ]] && die "Missing R2 for pair: $group"
    done

    declare -a _pids=() _ids=()
    for group in "${!pe_r1[@]}"; do
        quant_pe "$group" \
            "${pe_r1[$group]}" "${pe_r2[$group]}" \
            "${pe_acc1[$group]}" "${pe_acc2[$group]}" \
            "$TSV_FILE" &
        _pids+=($!) _ids+=("$group")
        if (( ${#_pids[@]} >= MAX_JOBS )); then
            wait "${_pids[0]}" || die "PE quant failed: ${_ids[0]}"
            _pids=("${_pids[@]:1}") _ids=("${_ids[@]:1}")
        fi
    done
    for i in "${!_pids[@]}"; do wait "${_pids[$i]}" || die "PE quant failed: ${_ids[$i]}"; done

    rm -f "${TSV_FILE}.lock"
    success "$TSV_BASE done | log: $LOG"

    unset se_order se_acc se_trimmed pe_r1 pe_r2 pe_acc1 pe_acc2
done

success "RNA-seq quantification pipeline complete!"