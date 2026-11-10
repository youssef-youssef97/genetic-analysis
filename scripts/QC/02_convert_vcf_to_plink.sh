
#!/bin/bash
set -euo pipefail

cd /home/youssef/Project/genetic-analysis

PLINK2="/home/youssef/Project/plink2_install/plink2"

VCF_DIR="data/raw/chromosomes/vcf"
OUTPUT_DIR="data/raw/chromosomes/plink"

PSAM="$VCF_DIR/1000G_samples.psam"
PSAM_MALES="$VCF_DIR/1000G_males.psam"

mkdir -p "$OUTPUT_DIR"

echo "======================================"
echo "1000 Genomes - VCF to PLINK 2"
echo "======================================"

convert_chromosome() {

    local chr="$1"
    local vcf="$2"
    local output="$OUTPUT_DIR/chr${chr}"
    local sample_file="$PSAM"

    if [[ ! -f "$vcf" ]]; then
        echo "ERROR: Missing VCF for chr${chr}"
        exit 1
    fi

    if [[ "$chr" == "Y" ]]; then
        sample_file="$PSAM_MALES"
    fi

    if [[ ! -f "$sample_file" ]]; then
        echo "ERROR: Missing sample file: $sample_file"
        exit 1
    fi

    echo ""
    echo "Converting chromosome ${chr}..."

    if [[ "$chr" == "X" ]]; then
        "$PLINK2" \
            --vcf "$vcf" \
            --psam "$sample_file" \
            --double-id \
            --split-par b37 \
            --set-all-var-ids '@:#:$r:$a' \
            --var-id-multi '@:#:$r' \
            --new-id-max-allele-len 700 \
            --make-pgen \
            --out "$output"

    else
        "$PLINK2" \
            --vcf "$vcf" \
            --psam "$sample_file" \
            --double-id \
            --set-all-var-ids '@:#:$r:$a' \
            --var-id-multi '@:#:$r' \
            --new-id-max-allele-len 700 \
            --make-pgen \
            --out "$output"
    fi

    echo "chr${chr} completed."
}

# Autosomes
# for chr in {1..22}; do
#     VCF="$VCF_DIR/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
#     convert_chromosome "$chr" "$VCF"
# done

# Chromosome X
# convert_chromosome "X" \
#     "$VCF_DIR/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz"

# Chromosome Y
convert_chromosome "Y" \
    "$VCF_DIR/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz"

echo ""
echo "======================================"
echo "All 24 chromosomes converted successfully."
echo "Output: $OUTPUT_DIR"
echo "======================================"
