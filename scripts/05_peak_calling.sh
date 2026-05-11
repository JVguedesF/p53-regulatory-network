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
Q_VALUE="0.05"
GENOME_SIZE="hs"

_CPUS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
MAX_JOBS="$_CPUS"

call_peaks() {
    local sample_id="$1" chip_bam="$2" input_bam="$3" acc="$4" tsv="$5" peaks_out="$6"
    local narrowpeak="${peaks_out}/${sample_id}_peaks.narrowPeak"

    if [[ -s "$narrowpeak" ]]; then
        log "Skipping $sample_id — already processed."
        local n_peaks; n_peaks=$(wc -l < "$narrowpeak")
        tsv_update "$tsv" "$acc" "Peak_File=$narrowpeak" "N_Peaks=$n_peaks"
        return 0
    fi

    samtools quickcheck "$chip_bam"  || die "Corrupt BAM: $chip_bam"
    samtools quickcheck "$input_bam" || die "Corrupt BAM: $input_bam"
    [[ -f "${chip_bam}.bai" ]]       || die "Missing index: ${chip_bam}.bai"

    write_log_header "MACS2" "$sample_id" "$(basename "$chip_bam") vs $(basename "$input_bam")"

    local tmp; tmp=$(mktemp)
    macs2 callpeak \
        -t "$chip_bam" \
        -c "$input_bam" \
        -f BAM \
        -g "$GENOME_SIZE" \
        -n "$sample_id" \
        --call-summits \
        -q "$Q_VALUE" \
        --outdir "$peaks_out" > "$tmp" 2>&1 \
    || { echo "[$(date '+%H:%M:%S')] ERROR: $sample_id" >> "$LOG"; cat "$tmp" >> "$LOG"; rm -f "$tmp"; exit 1; }

    { flock 200
      cat "$tmp"
      echo "[$(date '+%H:%M:%S')] END: $sample_id"
    } >> "$LOG" 200>> "$LOG"
    rm -f "$tmp"

    if [[ ! -f "$narrowpeak" ]]; then
        warn "MACS2 produced no .narrowPeak for $sample_id"
        return 0
    fi

    local n_peaks; n_peaks=$(wc -l < "$narrowpeak")
    [[ "$n_peaks" -eq 0 ]] && warn "0 peaks for $sample_id — inspect ${peaks_out}/${sample_id}_peaks.xls"

    tsv_update "$tsv" "$acc" "Peak_File=$narrowpeak" "N_Peaks=$n_peaks"
    success "$sample_id — $n_peaks peaks"
}

for cmd in macs2 samtools python3 flock; do
    command -v "$cmd" &>/dev/null || die "'$cmd' not found in PATH."
done
python3 -c "import src.tsv_updater" 2>/dev/null \
    || die "src/tsv_updater.py not importable — run from project root."

collect_tsv_files "$@"

# shellcheck disable=SC2153
for TSV_FILE in "${TSV_FILES[@]}"; do
    [[ -f "$TSV_FILE" ]] || { warn "$TSV_FILE not found — skipping."; continue; }
    TSV_BASE="${TSV_FILE##*/}"; TSV_BASE="${TSV_BASE%.tsv}"

    declare -A _bam_of_acc=()
    has_chipseq=false

    while IFS=$'\x01' read -r _cond acc sample_type _pair _pe _rep \
                                 _input_acc _fastq _dl _sz _md5 \
                                 _qc _trimmed bam_path _rest; do
        [[ -z "$acc" || -z "$sample_type" ]] && continue
        case "${sample_type,,}" in
            ip|input) has_chipseq=true; _bam_of_acc["$acc"]="$bam_path" ;;
            rna) continue ;;
        esac
    done < <(tail -n +2 "$TSV_FILE" | tr $'\t' $'\x01')

    if ! "$has_chipseq"; then
        log "No ChIP-seq rows in $TSV_FILE — skipping."
        continue
    fi

    PEAKS_BASE="${OUT_DIR}/chipseq/${TSV_BASE}/05_peaks"
    LOG="${LOG_DIR}/chipseq/${TSV_BASE}_peak_calling.log"
    mkdir -p "$PEAKS_BASE" "$(dirname "$LOG")"

    log "Processing: $TSV_FILE | jobs: $MAX_JOBS"

    declare -a _pids=() _ids=()

    while IFS=$'\x01' read -r cond acc sample_type _pair _pe rep \
                                 input_acc _fastq _dl _sz _md5 \
                                 _qc _trimmed chip_bam _rest; do
        [[ -z "$acc" || -z "$sample_type" ]] && continue

        case "${sample_type,,}" in
            rna|input) continue ;;
            ip) ;;
            *) die "Unknown Sample_Type='$sample_type' for $acc" ;;
        esac

        [[ -f "$chip_bam" ]] \
            || die "IP BAM not found for $acc: '$chip_bam' (run 04 first?)"
        [[ -n "$input_acc" ]] \
            || die "Input_Accession is empty for IP $acc"

        local_input="${_bam_of_acc[$input_acc]:-}"
        [[ -n "$local_input" ]] || die "Accession '$input_acc' not found in TSV"
        [[ -f "$local_input" ]] || die "Input BAM not found: '$local_input'"

        sample_id="${cond}_Rep${rep}_${acc}"
        sample_out="${PEAKS_BASE}/${sample_id}"
        mkdir -p "$sample_out"

        call_peaks "$sample_id" "$chip_bam" "$local_input" "$acc" "$TSV_FILE" "$sample_out" &
        _pids+=("$!") _ids+=("$sample_id")

        if (( ${#_pids[@]} >= MAX_JOBS )); then
            wait "${_pids[0]}" || die "Peak calling failed: ${_ids[0]}"
            _pids=("${_pids[@]:1}") _ids=("${_ids[@]:1}")
        fi
    done < <(tail -n +2 "$TSV_FILE" | tr $'\t' $'\x01')

    for i in "${!_pids[@]}"; do
        wait "${_pids[$i]}" || die "Peak calling failed: ${_ids[$i]}"
    done

    rm -f "${TSV_FILE}.lock"
    success "$TSV_BASE done | log: $LOG"

    unset _bam_of_acc
done

success "Peak calling pipeline complete!"