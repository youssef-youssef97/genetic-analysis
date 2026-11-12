PLINK2="/home/youssef/Project/plink2_install/plink2"

INPUT="snp_array/data/processed/matched/1kg_omni_matched"

OUT="snp_array/results/PCA/1kg_omni_matched"

"$PLINK2" \
    --bfile "$INPUT" \
    --extract "${OUT}.prune.in" \
    --pca 20 \
    --out "$OUT"