#!/bin/bash
set -euo pipefail
shopt -s nullglob

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$(dirname "$SCRIPT_DIR")"

(( BASH_VERSINFO[0] >= 4 )) || { echo "Bash 4+ required (found $BASH_VERSION)" >&2; exit 1; }

TSV_DIR="data/tsv"
OUT_DIR="pipeline_outputs"
Q_VALUE="0.05"
GENOME_SIZE="hs"
MAX_JOBS=$(nproc 2>/dev/null || echo 4)

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; NC='\033[0m'

log()     { echo "[$(date '+%H:%M:%S')] $*"; [[ -n "${LOG:-}" ]] && echo "[$(date '+%H:%M:%S')] $*" >> "$LOG"; }
die()     { echo -e "${RED}FATAL ERROR:${NC} $*" >&2; [[ -n "${LOG:-}" ]] && echo "FATAL ERROR: $*" >> "$LOG"; exit 1; }
warn()    { echo -e "${YELLOW}WARNING:${NC} $*" >&2; [[ -n "${LOG:-}" ]] && echo "WARNING: $*" >> "$LOG"; }
success() { echo -e "${GREEN}OK${NC} $*"; [[ -n "${LOG:-}" ]] && echo "OK $*" >> "$LOG"; }

tsv_update() {
    local tsv="$1" acc="$2"; shift 2
    flock -x "${tsv}.lock" python src/tsv_updater.py "$tsv" "$acc" "$@"
}

write_log_header() {
    local sample_id="$1"
    { flock 200
      printf '\n%s\n' '====================================================================='
      echo "[$(date '+%H:%M:%S')] START: $sample_id"
      printf '%s\n' '---------------------------------------------------------------------'
    } >> "$LOG" 200>> "$LOG"
}

call_peaks() {
    local sample_id="$1" chip_bam="$2" input_bam="$3" acc="$4" tsv="$5" peaks_out="$6"

    local narrowpeak="${peaks_out}/${sample_id}_peaks.narrowPeak"

    if [[ -s "$narrowpeak" ]]; then
        log "Skipping $sample_id — already processed."
        return 0
    fi

    samtools quickcheck "$chip_bam"  || die "Corrupt BAM: $chip_bam"
    samtools quickcheck "$input_bam" || die "Corrupt BAM: $input_bam"
    [[ -f "${chip_bam}.bai" ]]       || die "Missing index: ${chip_bam}.bai"

    write_log_header "$sample_id"

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

    local n_peaks
    n_peaks=$(wc -l < "$narrowpeak")

    if [[ "$n_peaks" -eq 0 ]]; then
        warn "0 peaks for $sample_id — inspect ${peaks_out}/${sample_id}_peaks.xls"
    fi

    tsv_update "$tsv" "$acc" "Peak_File=$narrowpeak" "N_Peaks=$n_peaks"
    success "$sample_id — $n_peaks peaks"
}

for cmd in macs2 samtools python; do
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

    PEAKS_BASE="${OUT_DIR}/chipseq/${TSV_BASE}/05_peaks"
    LOG="logs/chipseq/${TSV_BASE}_peak_calling.log"
    mkdir -p "$PEAKS_BASE" "logs/chipseq"

    log "Processing: $TSV_FILE"

    declare -A bam_of_acc=()

    while IFS=$'\x01' read -r _cond acc sample_type _pair _pe _rep \
                                 _input_acc _fastq _dl _sz _md5 \
                                 _qc _trimmed bam_path _rest; do
        [[ -z "$acc" || -z "$sample_type" ]] && continue
        case "${sample_type,,}" in
            ip|input) bam_of_acc["$acc"]="$bam_path" ;;
        esac
    done < <(tail -n +2 "$TSV_FILE" | tr $'\t' $'\x01')

    declare -a _pids=() _ids=()

    while IFS=$'\x01' read -r cond acc sample_type _pair _pe rep \
                                 input_acc _fastq _dl _sz _md5 \
                                 _qc _trimmed chip_bam _rest; do
        [[ -z "$acc" || -z "$sample_type" ]] && continue

        case "${sample_type,,}" in
            rna)   log "Skipping RNA row: $acc"; continue ;;
            input) continue ;;
            ip)    ;;
            *)     die "Unknown Sample_Type='$sample_type' for $acc" ;;
        esac

        [[ -f "$chip_bam" ]] \
            || die "BAM not found for $acc: '$chip_bam' (did script 03 run?)"
        [[ -n "$input_acc" ]] \
            || die "Input_Accession is empty for IP $acc"

        local_input="${bam_of_acc[$input_acc]:-}"
        [[ -n "$local_input" ]] \
            || die "Accession '$input_acc' not found in TSV"
        [[ -f "$local_input" ]] \
            || die "Input BAM not found: '$local_input'"

        sample_id="${cond}_Rep${rep}_${acc}"
        sample_out="${PEAKS_BASE}/${sample_id}"
        mkdir -p "$sample_out"

        call_peaks "$sample_id" "$chip_bam" "$local_input" "$acc" "$TSV_FILE" "$sample_out" &
        _pids+=($!) _ids+=("$sample_id")

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

    unset bam_of_acc
done

success "Peak calling pipeline complete!"