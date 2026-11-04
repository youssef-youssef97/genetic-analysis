#!/bin/bash

set -euo pipefail

PLINK="/home/youssef/Project/plink"

INPUT="data/raw/chr14/plink/plink"
OUTPUT="results/QC/chr14_initial"

mkdir -p results/QC

echo "Starting initial QC for chromosome 14..."

"$PLINK" \
    --bfile "$INPUT" \
    --missing \
    --freq \
    --hardy \
    --out "$OUTPUT"

echo "Initial QC completed successfully."