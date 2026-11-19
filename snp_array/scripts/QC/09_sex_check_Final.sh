PLINK2="/home/youssef/Project/plink2_install/plink2"

INPUT="snp_array/data/processed/sex_check/omni_XY_matched"
OUT="snp_array/results/QC/SexCheck/sex_check_thresholds"

mkdir -p "$(dirname "$OUT")"

"$PLINK2" \
    --bfile "$INPUT" \
    --check-sex \
    max-female-xf=0.60 \
    min-male-xf=0.80 \
    max-female-yrate=0.30 \
    min-male-yrate=0.45 \
    --out "$OUT"