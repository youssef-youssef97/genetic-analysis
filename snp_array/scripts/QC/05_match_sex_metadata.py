
import pandas as pd
from pathlib import Path

BASE = Path("snp_array")
DATA = BASE / "data/raw/plink"
MATCHED = BASE / "data/processed/matched"
OUT = BASE / "data/processed/sex_check"

OUT.mkdir(parents=True, exist_ok=True)

# Read XY FAM
fam = pd.read_csv(
    DATA / "omni_XY.fam",
    sep=r"\s+",
    header=None,
    names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"]
)

# Read matched metadata
metadata = pd.read_csv(
    MATCHED / "matched_metadata.psam",
    sep=r"\s+"
)

# Match samples using IID
merged = fam[["FID", "IID"]].merge(
    metadata[["IID", "SEX"]],
    on="IID",
    how="inner",
    validate="one_to_one"
)

# Check matching
print("XY samples:", len(fam))
print("Matched samples:", len(merged))
print("Unmatched XY samples:", len(fam) - len(merged))

# Create PLINK update-sex file: FID IID SEX
sex_file = merged[["FID", "IID", "SEX"]].copy()

sex_file["SEX"] = sex_file["SEX"].astype(int)

sex_file.to_csv(
    OUT / "update_sex.txt",
    sep="\t",
    index=False,
    header=False
)

print("\nSex distribution:")
print(sex_file["SEX"].value_counts())

print("\nSaved:", OUT / "update_sex.txt")
