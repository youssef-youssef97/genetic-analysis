
#!/bin/bash

set -euo pipefail

PLINK="/home/youssef/Project/plink"

INPUT="data/raw/chr14/plink/plink"
OUTPUT="data/processed/chr14/plink/chr14_qc"

mkdir -p data/processed/chr14/plink

"$PLINK" \
    --bfile "$INPUT" \
    --mind 0.01 \
    --geno 0.01 \
    --make-bed \
    --out "$OUTPUT"

echo "QC filtering completed successfully."
