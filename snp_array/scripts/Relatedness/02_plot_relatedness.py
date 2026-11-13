
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
    "Some samples have missing SuperPop metadata!"

assert samples["Population"].notna().all(), \
    "Some samples have missing Population metadata!"

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
# 4. Calculate Superpopulation boundaries
# --------------------------------------------------

superpop_groups = samples.groupby(
    "SuperPop", sort=False
).indices

superpop_boundaries = []
superpop_centers = []
superpop_labels = []

start = 0

for group, indices in superpop_groups.items():
    end = start + len(indices)

    superpop_boundaries.append(end)
    superpop_centers.append((start + end - 1) / 2)
    superpop_labels.append(group)

    start = end

# --------------------------------------------------
# 5. Calculate Population boundaries and centers
# --------------------------------------------------

population_groups = samples.groupby(
    ["SuperPop", "Population"],
    sort=False
).indices

population_boundaries = []
population_centers = []
population_labels = []

start = 0

for (superpop, population), indices in population_groups.items():
    end = start + len(indices)

    population_boundaries.append(end)
    population_centers.append((start + end - 1) / 2)
    population_labels.append(population)

    start = end

# --------------------------------------------------
# Helper function for heatmaps
# --------------------------------------------------

def plot_heatmap(
    data,
    filename,
    title,
    vmin,
    vmax,
    cmap_name="coolwarm",
    missing_color="white"
):

    fig, ax = plt.subplots(figsize=(22, 20))

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

    # Population boundaries: thin gray lines
    for boundary in population_boundaries[:-1]:
        ax.axhline(
            boundary - 0.5,
            color="gray",
            linewidth=0.45,
            alpha=0.8
        )
        ax.axvline(
            boundary - 0.5,
            color="gray",
            linewidth=0.45,
            alpha=0.8
        )

    # Superpopulation boundaries: thick black lines
    for boundary in superpop_boundaries[:-1]:
        ax.axhline(
            boundary - 0.5,
            color="black",
            linewidth=1.2
        )
        ax.axvline(
            boundary - 0.5,
            color="black",
            linewidth=1.2
        )

    # Population labels
    ax.set_xticks(population_centers)
    ax.set_yticks(population_centers)

    ax.set_xticklabels(
        population_labels,
        rotation=90,
        fontsize=7
    )

    ax.set_yticklabels(
        population_labels,
        fontsize=7
    )

    ax.set_xlabel(
        "Population",
        fontsize=11,
        labelpad=12
    )

    ax.set_ylabel(
        "Population",
        fontsize=11,
        labelpad=12
    )

    ax.set_title(
        title,
        fontsize=15,
        pad=18
    )

    # Superpopulation labels on top
    ax_top = ax.secondary_xaxis("top")
    ax_top.set_xticks(superpop_centers)
    ax_top.set_xticklabels(superpop_labels, fontsize=10)
    ax_top.set_xlabel("Superpopulation", labelpad=10)

    # Superpopulation labels on right
    ax_right = ax.secondary_yaxis("right")
    ax_right.set_yticks(superpop_centers)
    ax_right.set_yticklabels(superpop_labels, fontsize=10)
    ax_right.set_ylabel("Superpopulation", labelpad=12)

    # Colorbar
    cbar = fig.colorbar(
        im,
        ax=ax,
        fraction=0.035,
        pad=0.06
    )

    cbar.set_label("KING Kinship Coefficient")

    plt.tight_layout()

    output = RESULTS / filename

    plt.savefig(
        output,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    print("Saved:", output)


# --------------------------------------------------
# 6. Expanded Kinship Heatmap
# --------------------------------------------------

plot_heatmap(
    matrix,
    "full_kinship_heatmap_expanded.png",
    "Full KING Kinship Heatmap",
    vmin=-0.45,
    vmax=0.30
)

# --------------------------------------------------
# 7. Related Pairs Only Heatmap
# --------------------------------------------------

positive_matrix = matrix.copy()

positive_matrix[
    positive_matrix <= 0.0442
] = np.nan

plot_heatmap(
    positive_matrix,
    "related_pairs_heatmap.png",
    "KING Kinship Heatmap - Related Pairs Only",
    vmin=0.0442,
    vmax=0.30,
    cmap_name="magma",
    missing_color="#eeeeee"
)

# --------------------------------------------------
# 8. Read pairwise KING results
# --------------------------------------------------

pairs = pd.read_csv(
    RESULTS / "1kg_omni_matched.kin0",
    sep=r"\s+"
)

# --------------------------------------------------
# 9. Kinship Distribution Histogram
# --------------------------------------------------

plt.figure(figsize=(11, 6))

plt.hist(
    pairs["KINSHIP"],
    bins=100,
    edgecolor="black"
)

plt.axvline(
    0.0442,
    linestyle="--",
    label="Third-degree threshold"
)

plt.axvline(
    0.0884,
    linestyle="--",
    label="Second-degree threshold"
)

plt.axvline(
    0.177,
    linestyle="--",
    label="First-degree threshold"
)

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
# 10. Top 20 Related Pairs
# --------------------------------------------------

top = pairs.nlargest(20, "KINSHIP").copy()

top["Pair"] = top["IID1"] + " - " + top["IID2"]
top = top.sort_values("KINSHIP")

plt.figure(figsize=(11, 8))

plt.barh(
    top["Pair"],
    top["KINSHIP"]
)

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

plt.savefig(
    output,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("Saved:", output)

print("\nAll plots generated successfully!")
