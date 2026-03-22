#!/bin/bash
set -euo pipefail
shopt -s nullglob

# ── Anchor to project root ────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

# ── Bash version check ────────────────────────────────────────────────────────
(( BASH_VERSINFO[0] >= 4 )) || { echo "Bash 4+ required (found $BASH_VERSION)" >&2; exit 1; }

# ── Configuration ─────────────────────────────────────────────────────────────
RAW_DIR="data/raw"
TSV_DIR="data/tsv"
OUT_DIR="pipeline_outputs"
LOG_DIR="logs"
TRIM_QUALITY=20
MAX_JOBS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

# ── Colors + Utility functions ────────────────────────────────────────────────
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; NC='\033[0m'
log()     { echo "[$(date '+%H:%M:%S')] $*"; }
die()     { echo -e "${RED}FATAL ERROR:${NC} $*" >&2; exit 1; }
warn()    { echo -e "${YELLOW}WARNING:${NC} $*" >&2; }
success() { echo -e "${GREEN}OK${NC} $*"; }

# ── Helpers ───────────────────────────────────────────────────────────────────
write_log_header() {
    local mode="$1" id="$2" desc="$3"
    {
        echo ""
        echo "====================================================================="
        echo "[$(date '+%H:%M:%S')] START ${mode}: ${id} | ${desc}"
        echo "---------------------------------------------------------------------"
    } >> "$LOG"
}
move_reports() {
    find "$TRIMMED_OUT" -maxdepth 1 -name "${1}*_trimming_report.txt" \
        -exec mv {} "$TRIM_REPORTS"/ \;
}

# ── Trimming functions ────────────────────────────────────────────────────────
trim_se() {
    local id="$1" acc="$2" file="$3"
    local expected="${TRIMMED_OUT}/${acc}_trimmed.fq.gz"
    if [[ -f "$expected" ]] && gzip -t "$expected" 2>/dev/null; then
        log "Skipping SE: $id — intact."; return 0
    fi
    [[ -f "$expected" ]] && warn "$expected corrupted. Reprocessing..."

    log "SE trimming: $id ($acc)"
    write_log_header "SE" "$id" "$acc.fastq.gz"

    local logtmp; logtmp=$(mktemp)
    trim_galore --quality "$TRIM_QUALITY" \
        --output_dir "$TRIMMED_OUT" "$file" > /dev/null 2> "$logtmp" \
    || { echo "[$(date '+%H:%M:%S')] ERROR in $id" >> "$LOG"; cat "$logtmp" >> "$LOG"; rm -f "$logtmp"; exit 1; }

    {
        cat "$TRIMMED_OUT/${acc}"*_trimming_report.txt 2>/dev/null || true
        echo "[$(date '+%H:%M:%S')] END SE: $id"
    } >> "$LOG"
    move_reports "$acc"
    rm -f "$logtmp"
}

trim_pe() {
    local id="$1" r1="$2" r2="$3"
    local acc1 acc2 out1 out2
    acc1=$(basename "$r1" .fastq.gz)
    acc2=$(basename "$r2" .fastq.gz)
    out1="${TRIMMED_OUT}/${acc1}_val_1.fq.gz"
    out2="${TRIMMED_OUT}/${acc2}_val_2.fq.gz"

    if [[ -f "$out1" && -f "$out2" ]] && gzip -t "$out1" 2>/dev/null && gzip -t "$out2" 2>/dev/null; then
        log "Skipping PE: $id — intact."; return 0
    fi
    [[ -f "$out1" || -f "$out2" ]] && warn "Corrupted file(s) for $id. Reprocessing..."

    log "PE trimming: $id"
    write_log_header "PE" "$id" "${acc1}.fastq.gz + ${acc2}.fastq.gz"

    local logtmp; logtmp=$(mktemp)
    trim_galore --paired --quality "$TRIM_QUALITY" \
        --output_dir "$TRIMMED_OUT" "$r1" "$r2" > /dev/null 2> "$logtmp" \
    || { echo "[$(date '+%H:%M:%S')] ERROR in $id" >> "$LOG"; cat "$logtmp" >> "$LOG"; rm -f "$logtmp"; exit 1; }

    {
        cat "$TRIMMED_OUT/${acc1}"*_trimming_report.txt 2>/dev/null || true
        cat "$TRIMMED_OUT/${acc2}"*_trimming_report.txt 2>/dev/null || true
        echo "[$(date '+%H:%M:%S')] END PE: $id"
    } >> "$LOG"
    move_reports "$acc1"
    move_reports "$acc2"
    rm -f "$logtmp"
}

# ── Dependency check ──────────────────────────────────────────────────────────
for cmd in fastqc trim_galore multiqc gzip awk mktemp; do
    command -v "$cmd" &>/dev/null || die "'$cmd' not found in PATH."
done

# ── TSV selection ─────────────────────────────────────────────────────────────
if [[ $# -gt 0 ]]; then
    TSV_FILES=("$@")
else
    mapfile -t TSV_FILES < <(find "$TSV_DIR" -maxdepth 1 -name "*.tsv" | sort)
fi
[[ ${#TSV_FILES[@]} -eq 0 ]] && die "No .tsv files found in $TSV_DIR"

# ── Main loop ─────────────────────────────────────────────────────────────────
for TSV_FILE in "${TSV_FILES[@]}"; do
    [[ -f "$TSV_FILE" ]] || { warn "$TSV_FILE not found. Skipping."; continue; }
    TSV_BASE="${TSV_FILE##*/}"; TSV_BASE="${TSV_BASE%.tsv}"
    log "Starting: $TSV_FILE"

    unset s_id s_acc s_file s_paired pe_r1 pe_r2 se_order se_file_map se_acc_map missing

    ANALYSIS=""
    declare -a s_id=() s_acc=() s_file=() s_paired=()
    row_type=""
    line_num=1
    while IFS=$'\t' read -r cond acc paired rep _rest; do
        (( line_num++ ))
        [[ -z "$cond" || -z "$acc" || -z "$paired" || -z "$rep" ]] && \
            die "$TSV_FILE line $line_num: insufficient columns."
        case "${cond,,}" in
            *chip-seq*|*chip_seq*|*chipseq*) row_type="chipseq" ;;
            *rna-seq*|*rna_seq*|*rnaseq*)    row_type="rnaseq"  ;;
            *) die "$TSV_FILE line $line_num: invalid type '$cond'." ;;
        esac
        [[ -z "$ANALYSIS" ]] && ANALYSIS="$row_type"
        [[ "$ANALYSIS" != "$row_type" ]] && die "$TSV_FILE: mixed types!"
        file_path="${RAW_DIR}/${acc}.fastq.gz"
        [[ -f "$file_path" ]] || die "$TSV_FILE line $line_num: not found -> $file_path"
        s_id+=("${cond}_Rep${rep}")
        s_acc+=("$acc")
        s_file+=("$file_path")
        paired_upper="${paired^^}"
        if [[ "$paired_upper" == "FALSE" || "$paired_upper" == "F" || "$paired_upper" == "0" || "$paired_upper" == "NO" || "$paired_upper" == "N" ]]; then
            s_paired+=("FALSE")
        else
            s_paired+=("TRUE")
        fi
    done < <(tail -n +2 "$TSV_FILE")

    [[ -z "$ANALYSIS" || ${#s_id[@]} -eq 0 ]] && { warn "Empty TSV: $TSV_FILE. Skipping."; continue; }
    log "Type: $ANALYSIS | ${#s_id[@]} samples | parallelism: $MAX_JOBS"

    # ── Paths ─────────────────────────────────────────────────────────────────
    BASE="${OUT_DIR}/${ANALYSIS}/${TSV_BASE}"
    QC_RAW="${BASE}/01_qc/fastqc_raw"
    QC_TRIMMED="${BASE}/01_qc/fastqc_trimmed"
    TRIM_REPORTS="${BASE}/01_qc/trimming_reports"
    MULTIQC_OUT="${BASE}/01_qc/multiqc"
    TRIMMED_OUT="${BASE}/02_trimmed"
    LOG="${LOG_DIR}/${ANALYSIS}/${TSV_BASE}.log"
    mkdir -p "$QC_RAW" "$QC_TRIMMED" "$TRIM_REPORTS" "$MULTIQC_OUT" "$TRIMMED_OUT" "$(dirname "$LOG")"

    # ── Raw FastQC ────────────────────────────────────────────────────────────
    declare -a missing=()
    for f in "${s_file[@]}"; do
        base=$(basename "$f" .fastq.gz)
        [[ -f "$QC_RAW/${base}_fastqc.html" ]] || missing+=("$f")
    done
    if [[ ${#missing[@]} -gt 0 ]]; then
        log "Raw FastQC (${#missing[@]} files)..."
        fastqc --quiet --threads "$MAX_JOBS" "${missing[@]}" --outdir "$QC_RAW" \
            || die "Raw FastQC failed."
    else
        log "Raw FastQC already exists. Skipping."
    fi

    # ── SE / PE separation ────────────────────────────────────────────────────
    declare -A pe_r1=() pe_r2=()
    declare -a se_order=()
    declare -A se_file_map=() se_acc_map=()
    for i in "${!s_id[@]}"; do
        id="${s_id[$i]}"
        if [[ "${s_paired[$i]}" == "FALSE" ]]; then
            se_order+=("$id")
            se_file_map["$id"]="${s_file[$i]}"
            se_acc_map["$id"]="${s_acc[$i]}"
        else
            [[ "${s_acc[$i]}" == *_1 ]] && pe_r1["$id"]="${s_file[$i]}"
            [[ "${s_acc[$i]}" == *_2 ]] && pe_r2["$id"]="${s_file[$i]}"
        fi
    done

    # ── Single-End trimming ───────────────────────────────────────────────────
    declare -a _se_pids=() _se_ids=()
    for id in "${se_order[@]}"; do
        trim_se "$id" "${se_acc_map[$id]}" "${se_file_map[$id]}" &
        _se_pids+=($!) _se_ids+=("$id")
        if (( ${#_se_pids[@]} >= MAX_JOBS )); then
            wait "${_se_pids[0]}" || die "SE job failed: ${_se_ids[0]}"
            _se_pids=("${_se_pids[@]:1}") _se_ids=("${_se_ids[@]:1}")
        fi
    done
    for i in "${!_se_pids[@]}"; do
        wait "${_se_pids[$i]}" || die "SE job failed: ${_se_ids[$i]}"
    done

    # ── Paired-End trimming ───────────────────────────────────────────────────
    for id in "${!pe_r1[@]}"; do [[ -z "${pe_r2[$id]:-}" ]] && warn "Incomplete pair — missing _2 for $id"; done
    for id in "${!pe_r2[@]}"; do [[ -z "${pe_r1[$id]:-}" ]] && warn "Incomplete pair — missing _1 for $id"; done

    declare -a _pe_pids=() _pe_ids=()
    for id in "${!pe_r1[@]}"; do
        [[ -z "${pe_r2[$id]:-}" ]] && continue
        trim_pe "$id" "${pe_r1[$id]}" "${pe_r2[$id]}" &
        _pe_pids+=($!) _pe_ids+=("$id")
        if (( ${#_pe_pids[@]} >= MAX_JOBS )); then
            wait "${_pe_pids[0]}" || die "PE job failed: ${_pe_ids[0]}"
            _pe_pids=("${_pe_pids[@]:1}") _pe_ids=("${_pe_ids[@]:1}")
        fi
    done
    for i in "${!_pe_pids[@]}"; do
        wait "${_pe_pids[$i]}" || die "PE job failed: ${_pe_ids[$i]}"
    done

    # ── Trimmed FastQC ────────────────────────────────────────────────────────
    log "Trimmed FastQC..."
    declare -a trimmed_files=("$TRIMMED_OUT"/*.fq.gz)
    declare -a missing_trimmed=()
    for f in "${trimmed_files[@]}"; do
        base=$(basename "$f" .fq.gz)
        [[ -f "$QC_TRIMMED/${base}_fastqc.html" ]] || missing_trimmed+=("$f")
    done
    if [[ ${#missing_trimmed[@]} -gt 0 ]]; then
        fastqc --quiet --threads "$MAX_JOBS" "${missing_trimmed[@]}" --outdir "$QC_TRIMMED" \
            || die "Trimmed FastQC failed."
    else
        log "Trimmed FastQC already exists. Skipping."
    fi

    # ── MultiQC ───────────────────────────────────────────────────────────────
    log "MultiQC: $TSV_BASE..."
    multiqc "$QC_RAW" "$QC_TRIMMED" "$TRIM_REPORTS" \
        -o "$MULTIQC_OUT" \
        --filename "multiqc_${TSV_BASE}" \
        --title "${TSV_BASE} - ${ANALYSIS}" \
        --force > /dev/null 2>&1

    success "$TSV_BASE done | log: $LOG"
done

success "Pipeline complete!"