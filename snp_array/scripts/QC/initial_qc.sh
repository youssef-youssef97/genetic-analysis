
#!/bin/bash
set -euo pipefail

# Project directory
cd /home/youssef/Project/genetic-analysis

# PLINK 2.0 executable
PLINK2="/home/youssef/Project/plink2_install/plink2"

# SNP Array input files (BED/BIM/FAM)
INPUT="snp_array/data/raw/plink/1kg_omni_auto"

# QC output
QC_OUTPUT="snp_array/results/QC/1kg_omni_auto"

# Check input files
for ext in bed bim fam; do
    if [[ ! -f "${INPUT}.${ext}" ]]; then
        echo "Error: Missing ${INPUT}.${ext}"
        exit 1
    fi
done

# Create output directory
mkdir -p "$(dirname "$QC_OUTPUT")"

echo "Starting initial SNP Array QC..."

"$PLINK2" \
    --bfile "$INPUT" \
    --missing \
    --freq \
    --hardy \
    --het \
    --out "$QC_OUTPUT"

echo "Initial QC completed successfully!"
