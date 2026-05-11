#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$(dirname "$SCRIPT_DIR")"

(( BASH_VERSINFO[0] >= 4 )) || { echo "Bash 4+ required (found $BASH_VERSION)" >&2; exit 1; }

# shellcheck disable=SC1091
source src/pipeline_common.sh

# shellcheck disable=SC2034
TSV_DIR="data/tsv"
OUT_DIR="pipeline_outputs"
LOG_DIR="logs"

_CPUS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
THREADS="$_CPUS"
MAX_JOBS=1
SORT_MEM="2G"

declare -a FORWARD_ARGS=()
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --threads) THREADS="$2"; shift 2 ;;
        --jobs)    MAX_JOBS="$2"; shift 2 ;;
        --mem)     SORT_MEM="$2"; shift 2 ;;
        *)         FORWARD_ARGS+=("$1"); shift ;;
    esac
done
set -- "${FORWARD_ARGS[@]+"${FORWARD_ARGS[@]}"}"

dedup() {
    local id="$1" bam_in="$2" accs="$3" tsv="$4"
    local safe_id="${id// /_}"
    local bam_out="${DEDUP_OUT}/${safe_id}_dedup.bam"
    local metrics="${DEDUP_OUT}/reports/${safe_id}_metrics.txt"

    if [[ -f "$bam_out" ]] && samtools quickcheck "$bam_out" 2>/dev/null; then
        log "Skipping dedup: $id — intact."
        local acc; for acc in $accs; do tsv_update "$tsv" "$acc" "BAM_Path=$bam_out"; done
        return 0
    fi

    samtools quickcheck "$bam_in" || die "Corrupt input BAM: $bam_in"
    write_log_header "DEDUP" "$id" "$(basename "$bam_in")"

    local tmp tmp_sort
    tmp=$(mktemp)
    tmp_sort=$(mktemp -d)

    (
        samtools sort -n -m "$SORT_MEM" -@ "$THREADS" -T "${tmp_sort}/n" "$bam_in" |
        samtools fixmate -m - - |
        samtools sort -m "$SORT_MEM" -@ "$THREADS" -T "${tmp_sort}/c" - |
        samtools markdup -@ "$THREADS" -s -f "$metrics" - "$bam_out" &&
        samtools index -@ "$THREADS" "$bam_out"
    ) > "$tmp" 2>&1 \
    || { echo "[$(date '+%H:%M:%S')] ERROR $id" >> "$LOG"; cat "$tmp" >> "$LOG"; rm -rf "$tmp" "$tmp_sort"; exit 1; }

    { flock 200
      cat "$tmp"
      cat "$metrics" 2>/dev/null || true
      echo "[$(date '+%H:%M:%S')] END DEDUP: $id"
    } >> "$LOG" 200>> "$LOG"
    rm -rf "$tmp" "$tmp_sort"

    local dup_rate
    dup_rate=$(awk '/^DUPLICATE PRIMARY TOTAL:/{dr=$4} END{print dr+0}' "$metrics" 2>/dev/null || echo "N/A")

    local acc; for acc in $accs; do tsv_update "$tsv" "$acc" "BAM_Path=$bam_out"; done
    success "$id — duplicates removed (dup reads: ${dup_rate})"
}

for cmd in samtools python3 flock; do
    command -v "$cmd" &>/dev/null || die "'$cmd' not found in PATH."
done
python3 -c "import src.tsv_updater" 2>/dev/null \
    || die "src/tsv_updater.py not importable — run from project root."

collect_tsv_files "$@"

# shellcheck disable=SC2153
for TSV_FILE in "${TSV_FILES[@]}"; do
    [[ -f "$TSV_FILE" ]] || { warn "$TSV_FILE not found — skipping."; continue; }
    TSV_BASE="${TSV_FILE##*/}"; TSV_BASE="${TSV_BASE%.tsv}"

    has_chipseq=false
    while IFS=$'\x01' read -r _cond _acc sample_type _rest; do
        [[ -z "$sample_type" ]] && continue
        case "${sample_type,,}" in
            ip|input) has_chipseq=true; break ;;
        esac
    done < <(tail -n +2 "$TSV_FILE" | tr $'\t' $'\x01')

    if ! "$has_chipseq"; then
        log "No ChIP-seq rows in $TSV_FILE — skipping."
        continue
    fi

    DEDUP_OUT="${OUT_DIR}/chipseq/${TSV_BASE}/04_deduplication"
    LOG="${LOG_DIR}/chipseq/${TSV_BASE}_deduplication.log"
    mkdir -p "$DEDUP_OUT/reports" "$(dirname "$LOG")"

    log "Processing: $TSV_FILE | jobs: $MAX_JOBS | threads: $THREADS | mem: $SORT_MEM"

    declare -a _group_order=()
    declare -A _group_bam=() _group_accs=()

    while IFS=$'\x01' read -r cond acc sample_type _pair _pe rep \
                                 _input_acc _fastq _dl _sz _md5 \
                                 _qc _trimmed bam_path _rest; do
        [[ -z "$acc" || -z "$sample_type" ]] && continue

        case "${sample_type,,}" in
            rna)      log "Skipping RNA row: $acc"; continue ;;
            ip|input) ;;
            *)        die "Unknown Sample_Type='$sample_type' for $acc" ;;
        esac

        [[ -f "$bam_path" ]] \
            || die "BAM not found for $acc: '$bam_path' (run 03a first?)"

        sample_id="${cond}_${sample_type}_Rep${rep}"

        if [[ -z "${_group_bam[$sample_id]:-}" ]]; then
            _group_order+=("$sample_id")
            _group_bam["$sample_id"]="$bam_path"
            _group_accs["$sample_id"]="$acc"
        else
            _group_accs["$sample_id"]+=" $acc"
        fi
    done < <(tail -n +2 "$TSV_FILE" | tr $'\t' $'\x01')

    declare -a _pids=() _ids=()

    for sample_id in "${_group_order[@]}"; do
        dedup "$sample_id" "${_group_bam[$sample_id]}" "${_group_accs[$sample_id]}" "$TSV_FILE" &
        _pids+=("$!") _ids+=("$sample_id")

        if (( ${#_pids[@]} >= MAX_JOBS )); then
            wait "${_pids[0]}" || die "Dedup failed: ${_ids[0]}"
            _pids=("${_pids[@]:1}") _ids=("${_ids[@]:1}")
        fi
    done

    for i in "${!_pids[@]}"; do
        wait "${_pids[$i]}" || die "Dedup failed: ${_ids[$i]}"
    done

    unset _group_order _group_bam _group_accs

    rm -f "${TSV_FILE}.lock"
    success "$TSV_BASE done | log: $LOG"
done

success "Deduplication pipeline complete!"