
import pandas as pd
from pathlib import Path

BASE = Path("snp_array")
RESULTS = BASE / "results/Relatedness"
DATA = BASE / "data/raw/plink"

pairs = pd.read_csv(
    RESULTS / "1kg_omni_matched.kin0",
    sep=r"\s+"
)

metadata = pd.read_csv(
    DATA / "metadata.psam",
    sep=r"\s+"
).rename(columns={"#IID": "IID"})

pop = metadata.set_index("IID")[["SuperPop", "Population"]]

pairs["SuperPop1"] = pairs["IID1"].map(pop["SuperPop"])
pairs["SuperPop2"] = pairs["IID2"].map(pop["SuperPop"])

pairs["Population1"] = pairs["IID1"].map(pop["Population"])
pairs["Population2"] = pairs["IID2"].map(pop["Population"])

# Keep only pairs above the third-degree threshold
related = pairs[pairs["KINSHIP"] > 0.0442].copy()

related["Same_SuperPop"] = (
    related["SuperPop1"] == related["SuperPop2"]
)

related["Same_Population"] = (
    related["Population1"] == related["Population2"]
)

print("\nTotal related pairs:", len(related))

print("\nSame Superpopulation:")
print(related["Same_SuperPop"].value_counts())

print("\nSame Population:")
print(related["Same_Population"].value_counts())

print("\nRelated pairs by Population:")
print(
    related.groupby("Population1")
    .size()
    .sort_values(ascending=False)
)

print("\nRelated pairs by Superpopulation:")
print(
    related.groupby("SuperPop1")
    .size()
    .sort_values(ascending=False)
)

# Cross-population pair counts
print("\nPairs by Superpopulation combination:")
print(
    related.groupby(
        ["SuperPop1", "SuperPop2"]
    ).size().sort_values(ascending=False)
)

related.to_csv(
    RESULTS / "related_pairs_with_population.csv",
    index=False
)