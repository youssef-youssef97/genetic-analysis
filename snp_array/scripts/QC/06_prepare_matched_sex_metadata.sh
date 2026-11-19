
PLINK2="/home/youssef/Project/plink2_install/plink2"

DATA="snp_array/data/raw/plink/omni_XY"
OUT="snp_array/data/processed/sex_check/omni_XY_matched"
KEEP="snp_array/data/processed/sex_check/update_sex.txt"

# Create a separate keep file
awk '{print $1, $2}' "$KEEP" \
    > snp_array/data/processed/sex_check/keep_samples.txt

# Keep matched samples and update reported sex
"$PLINK2" \
    --bfile "$DATA" \
    --keep snp_array/data/processed/sex_check/keep_samples.txt \
    --update-sex "$KEEP" \
    --make-bed \
    --out "$OUT"
