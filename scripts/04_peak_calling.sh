#!/bin/bash
set -euo pipefail

CHIP_BAM="pipeline_outputs/test_peaks/test_chip.bam"
INPUT_BAM="pipeline_outputs/test_peaks/test_input.bam"
OUT_DIR="pipeline_outputs/test_peaks"

mkdir -p "$OUT_DIR"

macs2 callpeak \
  -t "$CHIP_BAM" \
  -c "$INPUT_BAM" \
  -f BAM \
  -g hs \
  -n test_p53 \
  --call-summits \
  --nomodel \
  --extsize 200 \
  -q 0.1 \
  --outdir "$OUT_DIR"

echo "=== files generated ==="
ls -lh "$OUT_DIR"

echo "=== peak count ==="
wc -l "$OUT_DIR/test_p53_peaks.narrowPeak"

echo "=== first lines ==="
head -3 "$OUT_DIR/test_p53_peaks.narrowPeak"