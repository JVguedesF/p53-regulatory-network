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
TRIM_QUALITY=20

_CPUS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
MAX_JOBS=$(( _CPUS / 2 < 1 ? 1 : _CPUS / 2 ))
TRIM_CORES=2

move_reports() {
    find "$TRIMMED_OUT" -maxdepth 1 -name "${1}*_trimming_report.txt" \
        -exec mv -f {} "$TRIM_REPORTS"/ \;
}

trim_se() {
    local id="$1" acc="$2" fastq="$3" tsv="$4"
    local trimmed="${TRIMMED_OUT}/${acc}_trimmed.fq.gz"

    if [[ -f "$trimmed" ]] && gzip -t "$trimmed" 2>/dev/null; then
        log "Skipping SE: $id — intact."
        tsv_update "$tsv" "$acc" "QC_Status=PASS" "Trimmed_Path=$trimmed"
        return 0
    fi
    [[ -f "$trimmed" ]] && warn "$trimmed corrupt — reprocessing."

    log "SE trimming: $id"
    write_log_header "SE" "$id" "$(basename "$fastq")"

    local tmp; tmp=$(mktemp)
    trim_galore --quality "$TRIM_QUALITY" --cores "$TRIM_CORES" \
        --output_dir "$TRIMMED_OUT" "$fastq" > "$tmp" 2>&1 \
    || { echo "[$(date '+%H:%M:%S')] ERROR $id" >> "$LOG"; cat "$tmp" >> "$LOG"; rm -f "$tmp"; exit 1; }

    { flock 200
      cat "$tmp"
      cat "$TRIMMED_OUT/${acc}"*_trimming_report.txt 2>/dev/null || true
      echo "[$(date '+%H:%M:%S')] END SE: $id"
    } >> "$LOG" 200>> "$LOG"

    move_reports "$acc"
    rm -f "$tmp"
    tsv_update "$tsv" "$acc" "QC_Status=PASS" "Trimmed_Path=$trimmed"
}

trim_pe() {
    local id="$1" r1="$2" r2="$3" acc1="$4" acc2="$5" tsv="$6"
    local out1="${TRIMMED_OUT}/${acc1}_val_1.fq.gz"
    local out2="${TRIMMED_OUT}/${acc2}_val_2.fq.gz"

    if [[ -f "$out1" && -f "$out2" ]] \
        && gzip -t "$out1" 2>/dev/null && gzip -t "$out2" 2>/dev/null; then
        log "Skipping PE: $id — intact."
        tsv_update "$tsv" "$acc1" "QC_Status=PASS" "Trimmed_Path=$out1"
        tsv_update "$tsv" "$acc2" "QC_Status=PASS" "Trimmed_Path=$out2"
        return 0
    fi
    [[ -f "$out1" || -f "$out2" ]] && warn "Partial output for $id — reprocessing."

    log "PE trimming: $id"
    write_log_header "PE" "$id" "$(basename "$r1") + $(basename "$r2")"

    local tmp; tmp=$(mktemp)
    trim_galore --paired --quality "$TRIM_QUALITY" --cores "$TRIM_CORES" \
        --output_dir "$TRIMMED_OUT" "$r1" "$r2" > "$tmp" 2>&1 \
    || { echo "[$(date '+%H:%M:%S')] ERROR $id" >> "$LOG"; cat "$tmp" >> "$LOG"; rm -f "$tmp"; exit 1; }

    { flock 200
      cat "$tmp"
      cat "$TRIMMED_OUT/${acc1}"*_trimming_report.txt 2>/dev/null || true
      cat "$TRIMMED_OUT/${acc2}"*_trimming_report.txt 2>/dev/null || true
      echo "[$(date '+%H:%M:%S')] END PE: $id"
    } >> "$LOG" 200>> "$LOG"

    move_reports "$acc1"
    move_reports "$acc2"
    rm -f "$tmp"
    tsv_update "$tsv" "$acc1" "QC_Status=PASS" "Trimmed_Path=$out1"
    tsv_update "$tsv" "$acc2" "QC_Status=PASS" "Trimmed_Path=$out2"
}

for cmd in fastqc trim_galore multiqc gzip python3 flock; do
    command -v "$cmd" &>/dev/null || die "'$cmd' not found in PATH."
done
python3 -c "import src.tsv_updater" 2>/dev/null \
    || die "src/tsv_updater.py not importable — run from project root."

collect_tsv_files "$@"

# shellcheck disable=SC2153
for TSV_FILE in "${TSV_FILES[@]}"; do
    [[ -f "$TSV_FILE" ]] || { warn "$TSV_FILE not found — skipping."; continue; }
    TSV_BASE="${TSV_FILE##*/}"; TSV_BASE="${TSV_BASE%.tsv}"

    for TARGET_ANALYSIS in chipseq rnaseq; do
        declare -a s_raw_fastq=() se_order=()
        declare -A se_acc=() se_fastq=() pe_r1=() pe_r2=() pe_acc1=() pe_acc2=()
        ANALYSIS=""

        line_num=1
        while IFS=$'\x01' read -r cond acc sample_type pair_id paired_end rep \
                                   _input_acc fastq_path _rest; do
            (( line_num++ ))
            [[ -z "$cond" || -z "$acc" || -z "$sample_type" || -z "$fastq_path" ]] \
                && die "Line $line_num: missing required columns."

            case "${sample_type,,}" in
                ip|input) row_type="chipseq" ;;
                rna)      row_type="rnaseq"  ;;
                *) die "Line $line_num: unknown Sample_Type='$sample_type'." ;;
            esac

            [[ "$row_type" != "$TARGET_ANALYSIS" ]] && continue
            ANALYSIS="$row_type"

            [[ -f "$fastq_path" ]] || die "Line $line_num: FASTQ not found → $fastq_path"
            s_raw_fastq+=("$fastq_path")

            group="${cond}_${sample_type}_Rep${rep}"
            if [[ "${paired_end,,}" == "true" ]]; then
                case "$pair_id" in
                    1) pe_r1["$group"]="$fastq_path"; pe_acc1["$group"]="$acc" ;;
                    2) pe_r2["$group"]="$fastq_path"; pe_acc2["$group"]="$acc" ;;
                    *) die "Line $line_num: unexpected Pair_ID='$pair_id'." ;;
                esac
            else
                se_order+=("$group")
                se_acc["$group"]="$acc"
                se_fastq["$group"]="$fastq_path"
            fi
        done < <(tail -n +2 "$TSV_FILE" | tr $'\t' $'\x01')

        [[ -z "$ANALYSIS" ]] && continue

        BASE="${OUT_DIR}/${ANALYSIS}/${TSV_BASE}"
        QC_RAW="${BASE}/01_qc/fastqc_raw"
        QC_TRIMMED="${BASE}/01_qc/fastqc_trimmed"
        TRIM_REPORTS="${BASE}/01_qc/trimming_reports"
        MULTIQC_OUT="${BASE}/01_qc/multiqc"
        TRIMMED_OUT="${BASE}/02_trimmed"
        LOG="${LOG_DIR}/${ANALYSIS}/${TSV_BASE}_qc_trimming.log"
        mkdir -p "$QC_RAW" "$QC_TRIMMED" "$TRIM_REPORTS" "$MULTIQC_OUT" \
                 "$TRIMMED_OUT" "$(dirname "$LOG")"

        log "Processing: $TSV_FILE | Type: $ANALYSIS | jobs: $MAX_JOBS | cores/job: $TRIM_CORES"

        declare -a missing_raw=()
        for f in "${s_raw_fastq[@]}"; do
            base=$(basename "$f" .fastq.gz)
            [[ -f "$QC_RAW/${base}_fastqc.html" ]] || missing_raw+=("$f")
        done
        (( ${#missing_raw[@]} > 0 )) && {
            fastqc --quiet --threads "$_CPUS" "${missing_raw[@]}" \
                --outdir "$QC_RAW" 2>/dev/null || die "Raw FastQC failed."
        }

        declare -a _pids=() _ids=()
        for group in "${se_order[@]}"; do
            trim_se "$group" "${se_acc[$group]}" "${se_fastq[$group]}" "$TSV_FILE" &
            _pids+=("$!") _ids+=("$group")
            if (( ${#_pids[@]} >= MAX_JOBS )); then
                wait "${_pids[0]}" || die "SE trim failed: ${_ids[0]}"
                _pids=("${_pids[@]:1}") _ids=("${_ids[@]:1}")
            fi
        done
        for i in "${!_pids[@]}"; do wait "${_pids[$i]}" || die "SE trim failed: ${_ids[$i]}"; done

        for group in "${!pe_r1[@]}"; do
            [[ -z "${pe_r2[$group]:-}" ]] && die "Missing R2 for pair: $group"
        done

        declare -a _pids=() _ids=()
        for group in "${!pe_r1[@]}"; do
            trim_pe "$group" "${pe_r1[$group]}" "${pe_r2[$group]}" \
                    "${pe_acc1[$group]}" "${pe_acc2[$group]}" "$TSV_FILE" &
            _pids+=("$!") _ids+=("$group")
            if (( ${#_pids[@]} >= MAX_JOBS )); then
                wait "${_pids[0]}" || die "PE trim failed: ${_ids[0]}"
                _pids=("${_pids[@]:1}") _ids=("${_ids[@]:1}")
            fi
        done
        for i in "${!_pids[@]}"; do wait "${_pids[$i]}" || die "PE trim failed: ${_ids[$i]}"; done

        declare -a missing_trimmed=()
        for f in "$TRIMMED_OUT"/*.fq.gz; do
            base=$(basename "$f" .fq.gz)
            [[ -f "$QC_TRIMMED/${base}_fastqc.html" ]] || missing_trimmed+=("$f")
        done
        (( ${#missing_trimmed[@]} > 0 )) && {
            fastqc --quiet --threads "$_CPUS" "${missing_trimmed[@]}" \
                --outdir "$QC_TRIMMED" 2>/dev/null || die "Trimmed FastQC failed."
        }

        multiqc "$QC_RAW" "$QC_TRIMMED" "$TRIM_REPORTS" \
            -o "$MULTIQC_OUT" --filename "multiqc_${TSV_BASE}" \
            --title "${TSV_BASE} — ${ANALYSIS}" --force >/dev/null 2>&1 \
            || warn "MultiQC failed."

        success "$TSV_BASE ($ANALYSIS) done | log: $LOG"

        unset s_raw_fastq se_order se_acc se_fastq pe_r1 pe_r2 pe_acc1 pe_acc2
    done

    rm -f "${TSV_FILE}.lock"
done

success "QC + Trimming pipeline complete!"