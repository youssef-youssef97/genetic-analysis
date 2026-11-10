
#!/bin/bash
set -euo pipefail

cd /home/youssef/Project/genetic-analysis

PLINK2="/home/youssef/Project/plink2_install/plink2"

INPUT="data/raw/chromosomes/plink"
OUTPUT="data/raw/chromosomes/plink/autosomes/autosomes"
QC_OUTPUT="results/QC/genomewide_initial"

mkdir -p "$(dirname "$OUTPUT")" results/QC

echo "Checking chromosome files..."

for chr in {1..22}; do
    for ext in pgen pvar psam; do
        if [[ ! -f "$INPUT/chr${chr}.${ext}" ]]; then
            echo "Error: Missing chr${chr}.${ext}"
            exit 1
        fi
    done
done

echo "All chromosome files found."
echo "Starting autosomes merging..."

MERGE_LIST="results/QC/merge_list.txt"
: > "$MERGE_LIST"

for chr in {2..22}; do
    echo "$INPUT/chr$chr" >> "$MERGE_LIST"
done

"$PLINK2" \
    --pfile "$INPUT/chr1" \
    --pmerge-list "$MERGE_LIST" \
    --make-pgen \
    --out "$OUTPUT"

echo "Merging completed."
echo "Starting genome-wide initial QC..."

"$PLINK2" \
    --pfile "$OUTPUT" \
    --missing \
    --freq \
    --hardy \
    --het \
    --out "$QC_OUTPUT"

echo "Genome-wide initial QC completed successfully."
