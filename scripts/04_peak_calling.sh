#!/bin/bash
set -euo pipefail
shopt -s nullglob

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$(dirname "$SCRIPT_DIR")"

(( BASH_VERSINFO[0] >= 4 )) || { echo "Bash 4+ required (found $BASH_VERSION)" >&2; exit 1; }

TSV_FILE="${1:-}"
OUT_DIR="pipeline_outputs/test_peaks"
Q_VALUE="0.1"
GENOME_SIZE="hs"

[[ -n "$TSV_FILE" ]] || { echo "Usage: $0 <tsv_file>" >&2; exit 1; }
[[ -f "$TSV_FILE" ]] || { echo "TSV not found: $TSV_FILE" >&2; exit 1; }

declare -A bam_of_acc=()

while IFS=$'\x01' read -r _cond acc sample_type _pair _pe _rep \
                             _input_acc _fastq _dl _sz _md5 \
                             _qc _trimmed bam_path _rest; do
    [[ -z "$acc" || -z "$sample_type" ]] && continue
    case "${sample_type,,}" in
        ip|input) bam_of_acc["$acc"]="$bam_path" ;;
    esac
done < <(tail -n +2 "$TSV_FILE" | tr $'\t' $'\x01')

chip_bam=""
input_bam=""
sample_id=""

while IFS=$'\x01' read -r cond acc sample_type _pair _pe rep \
                             input_acc _fastq _dl _sz _md5 \
                             _qc _trimmed bam_path _rest; do
    [[ -z "$acc" || -z "$sample_type" ]] && continue
    [[ "${sample_type,,}" == "ip" ]] || continue

    [[ -f "$bam_path" ]] \
        || { echo "BAM not found for $acc: '$bam_path'" >&2; exit 1; }

    [[ -n "$input_acc" ]] \
        || { echo "Input_Accession is empty for IP $acc" >&2; exit 1; }

    resolved_input="${bam_of_acc[$input_acc]:-}"
    [[ -n "$resolved_input" ]] \
        || { echo "Accession '$input_acc' not found in TSV." >&2; exit 1; }
    [[ -f "$resolved_input" ]] \
        || { echo "Input BAM not found: '$resolved_input'" >&2; exit 1; }

    chip_bam="$bam_path"
    input_bam="$resolved_input"
    sample_id="${cond}_Rep${rep}_${acc}"
    break
done < <(tail -n +2 "$TSV_FILE" | tr $'\t' $'\x01')

[[ -n "$chip_bam" ]] || { echo "No IP rows found in $TSV_FILE." >&2; exit 1; }

echo "Sample : $sample_id"
echo "ChIP   : $chip_bam"
echo "Input  : $input_bam"

mkdir -p "$OUT_DIR"

macs2 callpeak \
    -t "$chip_bam" \
    -c "$input_bam" \
    -f BAM \
    -g "$GENOME_SIZE" \
    -n "$sample_id" \
    --call-summits \
    --nomodel \
    --extsize 200 \
    -q "$Q_VALUE" \
    --outdir "$OUT_DIR"

echo "=== files generated ==="
ls -lh "$OUT_DIR"

echo "=== peak count ==="
wc -l "$OUT_DIR/${sample_id}_peaks.narrowPeak"

echo "=== first lines ==="
head -3 "$OUT_DIR/${sample_id}_peaks.narrowPeak"