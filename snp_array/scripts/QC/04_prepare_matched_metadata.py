
import pandas as pd
from pathlib import Path

BASE = Path("snp_array/data/raw/plink")
MATCHED = Path("snp_array/data/processed/matched")

# Load files
metadata = pd.read_csv(
    BASE / "metadata.psam",
    sep=r"\s+"
)

matched = pd.read_csv(
    MATCHED / "matched_samples.txt",
    sep=r"\s+",
    header=None,
    names=["FID", "IID"]
)

# Match metadata
matched_metadata = metadata[
    metadata["#IID"].isin(matched["IID"])
].copy()

# Preserve the order of matched samples
matched_metadata = (
    matched[["IID"]]
    .merge(
        matched_metadata.rename(columns={"#IID": "IID"}),
        on="IID",
        how="left",
        validate="one_to_one"
    )
)

# Save
output = MATCHED / "matched_metadata.psam"

matched_metadata.to_csv(
    output,
    sep="\t",
    index=False
)

print("Matched metadata created successfully.")
print("Samples:", len(matched_metadata))
print("Missing SEX:", matched_metadata["SEX"].isna().sum())
print("Output:", output)
