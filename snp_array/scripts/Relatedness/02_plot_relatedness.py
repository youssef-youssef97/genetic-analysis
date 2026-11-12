
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

BASE = Path("snp_array")
RESULTS = BASE / "results/Relatedness"
DATA = BASE / "data/raw/plink"

# --------------------------------------------------
# 1. Read KING matrix and sample IDs
# --------------------------------------------------

matrix_file = RESULTS / "1kg_omni_matched_matrix.king"
ids_file = RESULTS / "1kg_omni_matched_matrix.king.id"

matrix = np.loadtxt(matrix_file)
sample_ids = pd.read_csv(ids_file, sep=r"\s+")

print("Matrix shape:", matrix.shape)
print("Number of sample IDs:", len(sample_ids))

assert matrix.shape == (len(sample_ids), len(sample_ids)), \
    "Matrix dimensions do not match sample IDs!"

# --------------------------------------------------
# 2. Read population metadata
# --------------------------------------------------

metadata = pd.read_csv(DATA / "metadata.psam", sep=r"\s+")
metadata = metadata.rename(columns={"#IID": "IID"})

samples = sample_ids.merge(
    metadata[["IID", "SuperPop", "Population"]],
    on="IID",
    how="left",
    validate="one_to_one"
)

assert samples["SuperPop"].notna().all(), \
    "Some samples have missing population metadata!"

# --------------------------------------------------
# 3. Order samples by SuperPop and Population
# --------------------------------------------------

samples["_original_index"] = np.arange(len(samples))

samples = samples.sort_values(
    ["SuperPop", "Population", "IID"]
).reset_index(drop=True)

order = samples["_original_index"].to_numpy()
matrix = matrix[np.ix_(order, order)]

# Hide diagonal
np.fill_diagonal(matrix, np.nan)

# --------------------------------------------------
# 4. Calculate population boundaries
# --------------------------------------------------

group_positions = samples.groupby(
    "SuperPop", sort=False
).indices

boundaries = []
centers = []
labels = []

start = 0

for group, indices in group_positions.items():
    end = start + len(indices)

    boundaries.append(end)
    centers.append((start + end - 1) / 2)
    labels.append(group)

    start = end

# --------------------------------------------------
# Helper function for heatmaps
# --------------------------------------------------

def plot_heatmap(data, filename, title, vmin, vmax,
                 cmap_name="coolwarm", missing_color="white"):

    fig, ax = plt.subplots(figsize=(19, 17))

    cmap = plt.get_cmap(cmap_name).copy()
    cmap.set_bad(missing_color)

    im = ax.imshow(
        data,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        interpolation="nearest",
        aspect="equal"
    )

    for boundary in boundaries[:-1]:
        ax.axhline(boundary - 0.5, color="black", linewidth=0.7)
        ax.axvline(boundary - 0.5, color="black", linewidth=0.7)

    ax.set_xticks(centers)
    ax.set_yticks(centers)

    ax.set_xticklabels(labels, fontsize=11)
    ax.set_yticklabels(labels, fontsize=11)

    ax.set_xlabel("Samples grouped by Superpopulation")
    ax.set_ylabel("Samples grouped by Superpopulation")
    ax.set_title(title, fontsize=15, pad=18)

    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.04)
    cbar.set_label("KING Kinship Coefficient")

    plt.tight_layout()

    output = RESULTS / filename
    plt.savefig(output, dpi=300, bbox_inches="tight")
    plt.close()

    print("Saved:", output)

# --------------------------------------------------
# 5. Expanded Kinship Heatmap
# --------------------------------------------------

plot_heatmap(
    matrix,
    "full_kinship_heatmap_expanded.png",
    "Full KING Kinship Heatmap",
    vmin=-0.45,
    vmax=0.30
)

# --------------------------------------------------
# 6. Related Pairs Only Heatmap
# --------------------------------------------------

positive_matrix = matrix.copy()
positive_matrix[positive_matrix <= 0.0442] = np.nan

plot_heatmap(
    positive_matrix,
    "related_pairs_heatmap.png",
    "KING Kinship Heatmap - Related Pairs Only",
    vmin=0.0442,
    vmax=0.30,
    cmap_name="hot",
    missing_color="#eeeeee"
)

# --------------------------------------------------
# 7. Read pairwise KING results
# --------------------------------------------------

pairs = pd.read_csv(
    RESULTS / "1kg_omni_matched.kin0",
    sep=r"\s+"
)

# --------------------------------------------------
# 8. Kinship Distribution Histogram
# --------------------------------------------------

plt.figure(figsize=(11, 6))

plt.hist(
    pairs["KINSHIP"],
    bins=100,
    edgecolor="black"
)

plt.axvline(0.0442, linestyle="--", label="0.0442")
plt.axvline(0.0884, linestyle="--", label="0.0884")
plt.axvline(0.177, linestyle="--", label="0.177")

plt.xlabel("KING Kinship Coefficient")
plt.ylabel("Number of Sample Pairs")
plt.title("Distribution of Pairwise KING Kinship")
plt.legend()
plt.grid(alpha=0.2)

plt.tight_layout()

output = RESULTS / "kinship_distribution.png"
plt.savefig(output, dpi=300)
plt.close()

print("Saved:", output)

# --------------------------------------------------
# 9. Top 20 Related Pairs
# --------------------------------------------------

top = pairs.nlargest(20, "KINSHIP").copy()

top["Pair"] = top["IID1"] + " - " + top["IID2"]
top = top.sort_values("KINSHIP")

plt.figure(figsize=(11, 8))

plt.barh(top["Pair"], top["KINSHIP"])

plt.axvline(
    0.177,
    linestyle="--",
    label="First-degree threshold"
)

plt.xlabel("KING Kinship")
plt.ylabel("Sample Pair")
plt.title("Top 20 Most Related Sample Pairs")
plt.legend()
plt.tight_layout()

output = RESULTS / "top_20_related_pairs.png"
plt.savefig(output, dpi=300, bbox_inches="tight")
plt.close()

print("Saved:", output)

print("\nAll plots generated successfully!")
