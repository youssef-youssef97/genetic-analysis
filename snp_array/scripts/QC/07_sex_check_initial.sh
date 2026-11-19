PLINK2="/home/youssef/Project/plink2_install/plink2"

INPUT="snp_array/data/processed/sex_check/omni_XY_matched"

OUT="snp_array/results/QC/SexCheck/sex_check"

mkdir -p "$(dirname "$OUT")"

"$PLINK2" \
    --bfile "$INPUT" \
    --check-sex \
    --out "$OUT"