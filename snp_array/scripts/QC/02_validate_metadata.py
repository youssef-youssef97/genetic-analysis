
import pandas as pd
from pathlib import Path

# ============================================================
# PATHS
# ============================================================

BASE = Path("snp_array/data/raw/plink")

FAM_FILE = BASE / "1kg_omni_auto.fam"
META_FILE = BASE / "metadata.psam"

OUTPUT = Path("snp_array/results/QC/initial")
OUTPUT.mkdir(parents=True, exist_ok=True)

# ============================================================
# LOAD FILES
# ============================================================

print("Loading sample information...")

fam = pd.read_csv(
    FAM_FILE,
    sep=r"\s+",
    header=None,
    names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"]
)

metadata = pd.read_csv(META_FILE, sep=r"\s+")

metadata.rename(columns={"#IID": "IID"}, inplace=True)

# ============================================================
# VALIDATE METADATA
# ============================================================

required = ["IID", "PAT", "MAT", "SEX", "SuperPop", "Population"]

missing_columns = [
    col for col in required if col not in metadata.columns
]

if missing_columns:
    raise ValueError(f"Missing metadata columns: {missing_columns}")

print("\n========== BASIC INFORMATION ==========")

print(f"FAM samples: {len(fam):,}")
print(f"Metadata samples: {len(metadata):,}")

# ============================================================
# DUPLICATE CHECK
# ============================================================

print("\n========== DUPLICATE CHECK ==========")

fam_duplicates = fam["IID"].duplicated().sum()
meta_duplicates = metadata["IID"].duplicated().sum()

print("Duplicate FAM IDs:", fam_duplicates)
print("Duplicate Metadata IDs:", meta_duplicates)

if fam_duplicates or meta_duplicates:
    raise ValueError("Duplicate sample IDs detected!")

# ============================================================
# MISSING VALUES
# ============================================================

print("\n========== MISSING VALUES ==========")

print(metadata[required].isna().sum())

# ============================================================
# SAMPLE MATCHING
# ============================================================

print("\n========== SAMPLE MATCHING ==========")

fam_ids = set(fam["IID"])
meta_ids = set(metadata["IID"])

matched = fam_ids & meta_ids
missing_metadata = fam_ids - meta_ids
extra_metadata = meta_ids - fam_ids

print(f"Matched samples: {len(matched):,}")
print(f"Missing from metadata: {len(missing_metadata):,}")
print(f"Extra metadata samples: {len(extra_metadata):,}")

if missing_metadata:
    print("\nSamples missing from metadata:")
    print(sorted(missing_metadata)[:20])

if extra_metadata:
    print("\nExtra metadata samples (first 20):")
    print(sorted(extra_metadata)[:20])

# ============================================================
# POPULATION SUMMARY
# ============================================================

print("\n========== POPULATION DISTRIBUTION ==========")

matched_metadata = metadata[
    metadata["IID"].isin(fam_ids)
]

print("\nSuperPop:")
print(matched_metadata["SuperPop"].value_counts())

print("\nPopulation:")
print(matched_metadata["Population"].value_counts())

# ============================================================
# INVESTIGATE MISSING SAMPLES
# ============================================================

print("\n========== MISSING SAMPLE DETAILS ==========")

missing_samples = fam[
    ~fam["IID"].isin(meta_ids)
]

print(missing_samples.to_string(index=False))

missing_samples.to_csv(
    OUTPUT / "missing_metadata_samples.csv",
    index=False
)

# ============================================================
# INVESTIGATE MISSING SEX
# ============================================================

print("\n========== MISSING SEX INFORMATION ==========")

missing_sex = metadata[metadata["SEX"].isna()]

print(missing_sex.to_string(index=False))

matched_missing_sex = missing_sex[
    missing_sex["IID"].isin(fam_ids)
]

print(
    f"\nMatched samples with missing SEX: "
    f"{len(matched_missing_sex)}"
)

# ============================================================
# SAVE MATCHED SAMPLE IDS
# ============================================================

MATCHED_DIR = Path("snp_array/data/processed/matched")
MATCHED_DIR.mkdir(parents=True, exist_ok=True)

matched_samples = fam[
    fam["IID"].isin(meta_ids)
][["FID", "IID"]]

matched_samples.to_csv(
    MATCHED_DIR / "matched_samples.txt",
    sep="\t",
    index=False,
    header=False
)

print("\nMatched sample IDs saved:")
print(MATCHED_DIR / "matched_samples.txt")
print(f"Total: {len(matched_samples):,}")

# ============================================================
# FINAL STATUS
# ============================================================

print("\n========== FINAL STATUS ==========")

if not missing_metadata and len(matched) == len(fam):
    print("SUCCESS: All FAM samples have matching metadata.")
else:
    print("WARNING: Some samples are missing from metadata.")

print("\nMetadata validation completed.")
print(f"Missing sample report: {OUTPUT / 'missing_metadata_samples.csv'}")
