#!/bin/bash
set -euo pipefail
shopt -s nullglob

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$(dirname "$SCRIPT_DIR")"

(( BASH_VERSINFO[0] >= 4 )) || { echo "Bash 4+ required (found $BASH_VERSION)" >&2; exit 1; }

TSV_DIR="data/tsv"
OUT_DIR="pipeline_outputs"
LOG_DIR="logs"
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
    local id="$1" bam="$2"
    { flock 200
      printf '\n%s\n' '====================================================================='
      echo "[$(date '+%H:%M:%S')] START DEDUP: $id | $(basename "$bam")"
      printf '%s\n' '---------------------------------------------------------------------'
    } >> "$LOG" 200>> "$LOG"
}

dedup() {
    local id="$1" bam_in="$2" acc="$3" tsv="$4"
    local safe_id="${id// /_}"
    local bam_out="${DEDUP_OUT}/${safe_id}_${acc}_dedup.bam"
    local metrics="${DEDUP_OUT}/reports/${safe_id}_${acc}_metrics.txt"

    if [[ -f "$bam_out" ]] && samtools quickcheck "$bam_out" 2>/dev/null; then
        log "Skipping dedup: $id — intact."
        tsv_update "$tsv" "$acc" "BAM_Path=$bam_out"
        return 0
    fi

    samtools quickcheck "$bam_in" || die "Corrupt input BAM: $bam_in"
    write_log_header "$id" "$bam_in"

    local tmp; tmp=$(mktemp)
    (
        samtools sort -n -@ "$THREADS" "$bam_in" |
        samtools fixmate -m - - |
        samtools sort -@ "$THREADS" - |
        samtools markdup -@ "$THREADS" -s -f "$metrics" - "$bam_out" &&
        samtools index -@ "$THREADS" "$bam_out"
    ) > "$tmp" 2>&1 \
    || { echo "[$(date '+%H:%M:%S')] ERROR $id" >> "$LOG"; cat "$tmp" >> "$LOG"; rm -f "$tmp"; exit 1; }

    { flock 200
      cat "$tmp"
      cat "$metrics" 2>/dev/null || true
      echo "[$(date '+%H:%M:%S')] END DEDUP: $id"
    } >> "$LOG" 200>> "$LOG"
    rm -f "$tmp"

    local dup_rate
    dup_rate=$(awk '/^DUPLICATE READS:/{dr=$3} /^ESTIMATED LIBRARY SIZE:/{next} END{print dr+0}' "$metrics" 2>/dev/null || echo "N/A")

    tsv_update "$tsv" "$acc" "BAM_Path=$bam_out"
    success "$id — duplicates removed (dup reads: ${dup_rate})"
}

for cmd in samtools python; do
    command -v "$cmd" &>/dev/null || die "'$cmd' not found in PATH."
done
python -c "import src.tsv_updater" 2>/dev/null \
    || die "src/tsv_updater.py not importable — run from project root."

if [[ $# -gt 0 ]]; then
    TSV_FILES=("$@")
else
    mapfile -t TSV_FILES < <(find "$TSV_DIR" -maxdepth 1 -name "*.tsv" | sort)
fi
[[ ${#TSV_FILES[@]} -eq 0 ]] && die "No .tsv files found in $TSV_DIR"

for TSV_FILE in "${TSV_FILES[@]}"; do
    [[ -f "$TSV_FILE" ]] || { warn "$TSV_FILE not found — skipping."; continue; }
    TSV_BASE="${TSV_FILE##*/}"; TSV_BASE="${TSV_BASE%.tsv}"
    log "Processing: $TSV_FILE"

    DEDUP_OUT="${OUT_DIR}/chipseq/${TSV_BASE}/04_deduplication"
    LOG="${LOG_DIR}/chipseq/${TSV_BASE}_deduplication.log"
    mkdir -p "$DEDUP_OUT/reports" "$(dirname "$LOG")"

    declare -a _pids=() _ids=()

    while IFS=$'\x01' read -r cond acc sample_type _pair _pe rep \
                                 _input_acc _fastq _dl _sz _md5 \
                                 _qc _trimmed bam_path _rest; do
        [[ -z "$acc" || -z "$sample_type" ]] && continue

        case "${sample_type,,}" in
            rna) log "Skipping RNA row: $acc"; continue ;;
            ip|input) ;;
            *) die "Unknown Sample_Type='$sample_type' for $acc" ;;
        esac

        [[ -f "$bam_path" ]] \
            || die "BAM not found for $acc: '$bam_path' (did script 03 run?)"

        sample_id="${cond}_${sample_type}_Rep${rep}"

        dedup "$sample_id" "$bam_path" "$acc" "$TSV_FILE" &
        _pids+=($!) _ids+=("$sample_id")

        if (( ${#_pids[@]} >= MAX_JOBS )); then
            wait "${_pids[0]}" || die "Dedup failed: ${_ids[0]}"
            _pids=("${_pids[@]:1}") _ids=("${_ids[@]:1}")
        fi
    done < <(tail -n +2 "$TSV_FILE" | tr $'\t' $'\x01')

    for i in "${!_pids[@]}"; do
        wait "${_pids[$i]}" || die "Dedup failed: ${_ids[$i]}"
    done

    rm -f "${TSV_FILE}.lock"
    success "$TSV_BASE done | log: $LOG"
done

success "Deduplication pipeline complete!"