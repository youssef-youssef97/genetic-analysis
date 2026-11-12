
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

BASE = Path("snp_array")
PCA_DIR = BASE / "results/PCA"
DATA_DIR = BASE / "data/raw/plink"

# Read PCA results
eigenvec = pd.read_csv(
    PCA_DIR / "1kg_omni_matched.eigenvec",
    sep=r"\s+"
)

# Read population metadata
metadata = pd.read_csv(
    DATA_DIR / "metadata.psam",
    sep=r"\s+"
)

metadata = metadata.rename(columns={"#IID": "IID"})

# Merge PCA coordinates with population information
df = eigenvec.merge(
    metadata[["IID", "SuperPop", "Population"]],
    on="IID",
    how="left",
    validate="one_to_one"
)

print("Total samples:", len(df))
print("Missing SuperPop:", df["SuperPop"].isna().sum())
print("\nSuperPop distribution:")
print(df["SuperPop"].value_counts())

# Plot PC1 vs PC2
plt.figure(figsize=(10, 8))

for group in sorted(df["SuperPop"].dropna().unique()):
    subset = df[df["SuperPop"] == group]

    plt.scatter(
        subset["PC1"],
        subset["PC2"],
        label=group,
        alpha=0.7,
        s=22
    )

plt.xlabel("PC1")
plt.ylabel("PC2")
plt.title("PCA of 1000 Genomes Samples by Superpopulation")
plt.legend(title="Superpopulation")
plt.grid(alpha=0.2)
plt.tight_layout()

output = PCA_DIR / "PCA_superpopulation.png"
plt.savefig(output, dpi=300)
plt.close()

print("\nPlot saved to:", output)
