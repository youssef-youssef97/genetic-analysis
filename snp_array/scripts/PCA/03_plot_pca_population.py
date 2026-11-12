
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

BASE = Path("snp_array")
PCA_DIR = BASE / "results/PCA"
DATA_DIR = BASE / "data/raw/plink"

eigenvec = pd.read_csv(
    PCA_DIR / "1kg_omni_matched.eigenvec",
    sep=r"\s+"
)

metadata = pd.read_csv(
    DATA_DIR / "metadata.psam",
    sep=r"\s+"
).rename(columns={"#IID": "IID"})

df = eigenvec.merge(
    metadata[["IID", "SuperPop", "Population"]],
    on="IID",
    how="left",
    validate="one_to_one"
)

print("Missing Population:", df["Population"].isna().sum())
print("Number of populations:", df["Population"].nunique())

plt.figure(figsize=(14, 10))

populations = sorted(df["Population"].dropna().unique())
colors = plt.cm.tab20(
    [i / max(len(populations) - 1, 1)
     for i in range(len(populations))]
)

for population, color in zip(populations, colors):
    subset = df[df["Population"] == population]

    plt.scatter(
        subset["PC1"],
        subset["PC2"],
        label=population,
        color=color,
        alpha=0.75,
        s=20
    )

plt.xlabel("PC1")
plt.ylabel("PC2")
plt.title("PCA of 1000 Genomes Samples by Population")
plt.legend(
    title="Population",
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    fontsize=8
)

plt.grid(alpha=0.2)
plt.tight_layout()

output = PCA_DIR / "PCA_population.png"
plt.savefig(output, dpi=300, bbox_inches="tight")
plt.close()

print("Saved:", output)
