#!/bin/bash
set -euo pipefail
shopt -s nullglob

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -f "${SCRIPT_DIR}/config/config.yaml" ]]; then
    cd "$SCRIPT_DIR"
else
    cd "$(dirname "$SCRIPT_DIR")"
fi

(( BASH_VERSINFO[0] >= 4 )) || { echo "Bash 4+ required (found $BASH_VERSION)" >&2; exit 1; }

TSV_DIR="data/tsv"
OUT_DIR="pipeline_outputs"
LOG_DIR="logs"
TARGET="hg38"
GENOME_FASTA="data/genome/fasta/${TARGET}.fa"
INDEX="data/genome/index/bowtie2/genome_index_${TARGET}"
THREADS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
MAX_JOBS=1

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; NC='\033[0m'

log()     { echo "[$(date '+%H:%M:%S')] $*"; if [[ -n "${LOG:-}" ]]; then echo "[$(date '+%H:%M:%S')] $*" >> "$LOG"; fi; }
die()     { echo -e "${RED}FATAL ERROR:${NC} $*" >&2; if [[ -n "${LOG:-}" ]]; then echo "FATAL ERROR: $*" >> "$LOG"; fi; exit 1; }
warn()    { echo -e "${YELLOW}WARNING:${NC} $*" >&2; if [[ -n "${LOG:-}" ]]; then echo "WARNING: $*" >> "$LOG"; fi; }
success() { echo -e "${GREEN}OK${NC} $*"; if [[ -n "${LOG:-}" ]]; then echo "OK $*" >> "$LOG"; fi; }

tsv_update() {
    local tsv="$1" acc="$2"; shift 2
    flock -x "${tsv}.lock" python src/tsv_updater.py "$tsv" "$acc" "$@"
}

write_log_header() {
    local mode="$1" id="$2" desc="$3"
    { flock 200
      printf '\n%s\n' '====================================================================='
      echo "[$(date '+%H:%M:%S')] START ${mode}: ${id} | ${desc}"
      printf '%s\n' '---------------------------------------------------------------------'
    } >>"$LOG" 200>>"$LOG"
}

align_se() {
    local id="$1" acc="$2" trimmed="$3" tsv="$4"
    local safe_id="${id// /_}"
    local bam="${ALIGN_OUT}/${safe_id}_${acc}.bam"
    local report="${ALIGN_REPORTS}/${safe_id}_${acc}_report.txt"

    if [[ -f "$bam" ]] && samtools quickcheck "$bam" 2>/dev/null; then
        log "Skipping SE: $id — intact."; return 0
    fi
    log "SE alignment: $id"; write_log_header "ALIGN_SE" "$id" "$(basename "$trimmed")"

    local tmp; tmp=$(mktemp)
    (
        bowtie2 -p "$THREADS" -x "$INDEX" -U "$trimmed" 2>"$report" |
        samtools view -bS - |
        samtools sort -@ "$THREADS" -o "$bam" &&
        samtools index "$bam"
    ) >"$tmp" 2>&1 \
    || { echo "[$(date '+%H:%M:%S')] ERROR $id" >>"$LOG"; cat "$tmp" >>"$LOG"; rm -f "$tmp"; exit 1; }

    { flock 200
      cat "$tmp"
      cat "$report" 2>/dev/null || true
      echo "[$(date '+%H:%M:%S')] END ALIGN_SE: $id"
    } >>"$LOG" 200>>"$LOG"
    rm -f "$tmp"

    local rate; rate=$(grep "overall alignment rate" "$report" | grep -oP '[0-9.]+(?=%)' || echo "N/A")
    tsv_update "$tsv" "$acc" "BAM_Path=$bam" "Alignment_Rate=${rate}%"
}

align_pe() {
    local id="$1" r1="$2" r2="$3" acc1="$4" acc2="$5" acc_base="$6" tsv="$7"
    local safe_id="${id// /_}"
    local bam="${ALIGN_OUT}/${safe_id}_${acc_base}.bam"
    local report="${ALIGN_REPORTS}/${safe_id}_${acc_base}_report.txt"

    if [[ -f "$bam" ]] && samtools quickcheck "$bam" 2>/dev/null; then
        log "Skipping PE: $id — intact."; return 0
    fi
    log "PE alignment: $id"; write_log_header "ALIGN_PE" "$id" "${acc_base} (paired)"

    local tmp; tmp=$(mktemp)
    (
        bowtie2 -p "$THREADS" -x "$INDEX" -1 "$r1" -2 "$r2" 2>"$report" |
        samtools view -bS - |
        samtools sort -@ "$THREADS" -o "$bam" &&
        samtools index "$bam"
    ) >"$tmp" 2>&1 \
    || { echo "[$(date '+%H:%M:%S')] ERROR $id" >>"$LOG"; cat "$tmp" >>"$LOG"; rm -f "$tmp"; exit 1; }

    { flock 200
      cat "$tmp"
      cat "$report" 2>/dev/null || true
      echo "[$(date '+%H:%M:%S')] END ALIGN_PE: $id"
    } >>"$LOG" 200>>"$LOG"
    rm -f "$tmp"

    local rate; rate=$(grep "overall alignment rate" "$report" | grep -oP '[0-9.]+(?=%)' || echo "N/A")
    tsv_update "$tsv" "$acc1" "BAM_Path=$bam" "Alignment_Rate=${rate}%"
    tsv_update "$tsv" "$acc2" "BAM_Path=$bam" "Alignment_Rate=${rate}%"
}

for cmd in bowtie2 samtools python; do
    command -v "$cmd" &>/dev/null || die "'$cmd' not found in PATH."
done
python -c "import src.tsv_updater" 2>/dev/null \
    || die "src/tsv_updater.py not importable — run from project root."

mkdir -p "$LOG_DIR"
if [[ ! -f "${INDEX}.1.bt2" ]]; then
    log "Building genome index for $TARGET..."
    mkdir -p "$(dirname "$INDEX")"
    bowtie2-build "$GENOME_FASTA" "$INDEX" 2>&1 | tee -a "$LOG_DIR/genome_index.log" \
        || die "Failed to build genome index."
else
    log "Genome index found — skipping build."
fi

if [[ $# -gt 0 ]]; then
    TSV_FILES=("$@")
else
    mapfile -t TSV_FILES < <(find "$TSV_DIR" -maxdepth 1 -name "*.tsv" | sort)
fi
[[ ${#TSV_FILES[@]} -eq 0 ]] && die "No .tsv files found in $TSV_DIR"

for TSV_FILE in "${TSV_FILES[@]}"; do
    [[ -f "$TSV_FILE" ]] || { log "Skipping missing: $TSV_FILE"; continue; }
    TSV_BASE="${TSV_FILE##*/}"; TSV_BASE="${TSV_BASE%.tsv}"
    log "Processing: $TSV_FILE"

    ANALYSIS=""
    declare -a se_order=()
    declare -A se_acc=() se_trimmed=() \
               pe_r1=() pe_r2=() pe_acc1=() pe_acc2=() pe_acc_base=()

    line_num=1
    while IFS=$'\x01' read -r cond acc sample_type pair_id paired_end rep \
                               _input_acc _fastq_path _dl _size _md5 \
                               _qc trimmed_path _rest; do
        (( line_num++ ))
        [[ -z "$cond" || -z "$acc" || -z "$sample_type" ]] \
            && die "Line $line_num: missing required columns."

        case "${sample_type,,}" in
            ip|input) row_type="chipseq" ;;
            rna)
                log "Skipping RNA row $acc — STAR not implemented in this script."
                continue ;;
            *) die "Line $line_num: unknown Sample_Type='$sample_type'." ;;
        esac
        [[ -z "$ANALYSIS" ]] && ANALYSIS="$row_type"

        [[ -f "$trimmed_path" ]] \
            || die "Line $line_num: Trimmed_Path not found → '$trimmed_path' (run 02 first?)"

        group="${cond}_Rep${rep}"
        if [[ "${paired_end,,}" == "true" ]]; then
            case "$pair_id" in
                1) pe_r1["$group"]="$trimmed_path"; pe_acc1["$group"]="$acc"
                   pe_acc_base["$group"]="${acc%_1}" ;;
                2) pe_r2["$group"]="$trimmed_path"; pe_acc2["$group"]="$acc" ;;
                *) die "Line $line_num: unexpected Pair_ID='$pair_id'." ;;
            esac
        else
            se_order+=("$group")
            se_acc["$group"]="$acc"
            se_trimmed["$group"]="$trimmed_path"
        fi
    done < <(tail -n +2 "$TSV_FILE" | tr $'\t' $'\x01')

    [[ -z "$ANALYSIS" ]] && { log "No ChIP-seq rows in $TSV_FILE — skipping."; continue; }
    log "Type: $ANALYSIS | jobs: $MAX_JOBS"

    BASE="${OUT_DIR}/${ANALYSIS}/${TSV_BASE}"
    ALIGN_OUT="${BASE}/03_alignment"
    ALIGN_REPORTS="${ALIGN_OUT}/reports"
    LOG="${LOG_DIR}/${ANALYSIS}/${TSV_BASE}_alignment.log"
    mkdir -p "$ALIGN_OUT" "$ALIGN_REPORTS" "$(dirname "$LOG")"

    declare -a _pids=() _ids=()
    for group in "${se_order[@]}"; do
        align_se "$group" "${se_acc[$group]}" "${se_trimmed[$group]}" "$TSV_FILE" &
        _pids+=($!) _ids+=("$group")
        if (( ${#_pids[@]} >= MAX_JOBS )); then
            wait "${_pids[0]}" || die "SE alignment failed: ${_ids[0]}"
            _pids=("${_pids[@]:1}") _ids=("${_ids[@]:1}")
        fi
    done
    for i in "${!_pids[@]}"; do wait "${_pids[$i]}" || die "SE alignment failed: ${_ids[$i]}"; done

    for group in "${!pe_r1[@]}"; do
        [[ -z "${pe_r2[$group]:-}" ]] && die "Missing R2 for pair: $group"
    done

    declare -a _pids=() _ids=()
    for group in "${!pe_r1[@]}"; do
        align_pe "$group" \
            "${pe_r1[$group]}" "${pe_r2[$group]}" \
            "${pe_acc1[$group]}" "${pe_acc2[$group]}" \
            "${pe_acc_base[$group]}" "$TSV_FILE" &
        _pids+=($!) _ids+=("$group")
        if (( ${#_pids[@]} >= MAX_JOBS )); then
            wait "${_pids[0]}" || die "PE alignment failed: ${_ids[0]}"
            _pids=("${_pids[@]:1}") _ids=("${_ids[@]:1}")
        fi
    done
    for i in "${!_pids[@]}"; do wait "${_pids[$i]}" || die "PE alignment failed: ${_ids[$i]}"; done

    rm -f "${TSV_FILE}.lock"
    success "$TSV_BASE done | log: $LOG"

    unset se_order se_acc se_trimmed pe_r1 pe_r2 pe_acc1 pe_acc2 pe_acc_base
done

success "Alignment pipeline complete!"