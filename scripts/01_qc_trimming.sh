#!/bin/bash

set -e
set -o pipefail

mkdir -p data/processed results/qc/fastqc_raw results/qc/fastqc_trimmed logs results/reports/trimming

fastqc --quiet data/raw/*.fastq.gz --outdir results/qc/fastqc_raw

declare -A pair1
declare -A pair2

while IFS=$'\t' read -r cond acc paired rep rest; do
    if [[ "$cond" == "Condition" ]]; then continue; fi
    
    id="${cond}_Rep${rep}"
    file_path="data/raw/${acc}.fastq.gz"
    
    if [[ "$paired" == "False" || "$paired" == "FALSE" ]]; then
        echo "Starting Single-End trimming: $id ($acc)"
        trim_galore --quality 20 --fastqc --fastqc_args "--outdir results/qc/fastqc_trimmed/" --output_dir data/processed/ "$file_path" >> logs/trimming.log 2>&1 || { echo "CRITICAL ERROR: Trimming failed for sample $id ($acc)"; exit 1; }
        mv data/processed/*_trimming_report.txt results/reports/trimming/ 2>/dev/null || true
    else
        if [[ "$acc" == *_1 ]]; then
            pair1["$id"]="$file_path"
        elif [[ "$acc" == *_2 ]]; then
            pair2["$id"]="$file_path"
        fi
        
        if [[ -n "${pair1[$id]}" && -n "${pair2[$id]}" ]]; then
            echo "Starting Paired-End trimming: $id"
            trim_galore --paired --quality 20 --fastqc --fastqc_args "--outdir results/qc/fastqc_trimmed/" --output_dir data/processed/ "${pair1[$id]}" "${pair2[$id]}" >> logs/trimming.log 2>&1 || { echo "CRITICAL ERROR: Trimming failed for sample $id"; exit 1; }
            mv data/processed/*_trimming_report.txt results/reports/trimming/ 2>/dev/null || true
        fi
    fi
done < <(cat data/tsv/*.tsv)

echo "Processing complete. Check logs/trimming.log for details."