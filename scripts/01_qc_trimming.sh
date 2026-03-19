#!/bin/bash

shopt -s extglob nullglob

mkdir -p data/processed
mkdir -p results/qc/fastqc_raw
mkdir -p results/qc/fastqc_trimmed/

echo "Starting initial FastQC on raw data..."
fastqc data/raw/*.fastq.gz --outdir results/qc/fastqc_raw

mapfile -t samples < <(find data/raw -maxdepth 1 -name "*.fastq.gz" | sed -E 's/(_R[0-9]+|_[0-9]+)?\.fastq\.gz//' | sort | uniq)

for prefix in "${samples[@]}"
do
    files=( "$prefix"*fastq.gz )
    num_files=${#files[@]}
    base_name=$(basename "$prefix")

    if [ "$num_files" -eq 1 ]; then
        echo "--- Processing Single-End: $base_name ---"
        trim_galore --quality 20 --fastqc --fastqc_args "--outdir results/qc/fastqc_trimmed/" --output_dir data/processed/ "${files[0]}"

    elif [ "$num_files" -eq 2 ]; then
        echo "--- Processing Paired-End: $base_name ---"
        trim_galore --paired --quality 20 --fastqc --fastqc_args "--outdir results/qc/fastqc_trimmed/" --output_dir data/processed/ "${files[0]}" "${files[1]}"

    else
        echo "--- ALERT: $base_name has $num_files files. Check manually. ---"
    fi
done

echo "QC and Trimming process completed successfully!"