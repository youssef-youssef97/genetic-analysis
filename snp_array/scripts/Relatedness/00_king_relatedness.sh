PLINK2="/home/youssef/Project/plink2_install/plink2"
INPUT="snp_array/data/processed/matched/1kg_omni_matched"
SNPS="snp_array/results/PCA/1kg_omni_matched.prune.in"
OUT="snp_array/results/Relatedness/1kg_omni_matched"

mkdir -p snp_array/results/Relatedness

"$PLINK2" \
    --bfile "$INPUT" \
    --extract "$SNPS" \
    --make-king-table \
    --out "$OUT"

"$PLINK2" \
    --bfile "$INPUT" \
    --extract "$SNPS" \
    --make-king square \
    --out "${OUT}_matrix"