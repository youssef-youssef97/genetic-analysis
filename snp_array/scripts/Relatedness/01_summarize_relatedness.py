
import pandas as pd
from pathlib import Path

BASE = Path("snp_array/results/Relatedness")
FILE = BASE / "1kg_omni_matched.kin0"

df = pd.read_csv(FILE, sep=r"\s+")

print("Total pairs:", len(df))
print("\nKinship statistics:")
print(df["KINSHIP"].describe(percentiles=[
    0.50, 0.90, 0.95, 0.99, 0.999
]))

print("\nPairs above thresholds:")
for threshold in [0.0442, 0.0884, 0.177, 0.354]:
    count = (df["KINSHIP"] > threshold).sum()
    print(f"Kinship > {threshold}: {count}")

print("\nTop 20 pairs:")
top = df.nlargest(20, "KINSHIP")
print(
    top[["IID1", "IID2", "NSNP", "KINSHIP"]]
    .to_string(index=False)
)

top.to_csv(
    BASE / "top_related_pairs.csv",
    index=False
)
