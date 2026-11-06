
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# PATHS
# ============================================================

BASE = Path("results/QC/chr14_initial")
OUTPUT = Path("results/QC/chr14_qc_summary.txt")
PLOT_DIR = Path("results/QC/qc_plots")

PLOT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# LOAD DATA
# ============================================================

imiss = pd.read_csv(f"{BASE}.imiss", sep=r"\s+")
lmiss = pd.read_csv(f"{BASE}.lmiss", sep=r"\s+")
frq = pd.read_csv(f"{BASE}.frq", sep=r"\s+")
hwe = pd.read_csv(f"{BASE}.hwe", sep=r"\s+")
het = pd.read_csv(f"{BASE}.het", sep=r"\s+")

report = []

def add(text=""):
    report.append(str(text))


def summarize(name, x, thresholds=None):
    x = x.dropna()

    add(f"\n{'=' * 60}")
    add(name)
    add(f"Number of records: {len(x):,}")
    add(x.describe(
        percentiles=[0.01, 0.05, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99]
    ).to_string())

    if thresholds:
        for t in thresholds:
            count = (x < t).sum() if name == "MINOR ALLELE FREQUENCY" else (x > t).sum()
            add(f"Records at threshold {t:g}: {count:,} ({count / len(x) * 100:.4f}%)")


# ============================================================
# 1. MISSINGNESS STATISTICS
# ============================================================

summarize("SAMPLE MISSINGNESS", imiss["F_MISS"], [0.01, 0.02, 0.05])
summarize("VARIANT MISSINGNESS", lmiss["F_MISS"], [0.01, 0.02, 0.05])


# ============================================================
# 2. MAF STATISTICS
# ============================================================

maf = frq["MAF"]
summarize("MINOR ALLELE FREQUENCY", maf, [0.01, 0.05])


# ============================================================
# 3. HWE STATISTICS
# ============================================================

add(f"\n{'=' * 60}")
add("HARDY-WEINBERG EQUILIBRIUM")
add(f"Number of records: {len(hwe):,}")
add("\nTest categories:")
add(hwe["TEST"].value_counts().to_string())

for test in hwe["TEST"].unique():
    subset = hwe[hwe["TEST"] == test]
    add(f"\nTest: {test}")
    add(subset["P"].describe().to_string())

    for threshold in [1e-6, 1e-10]:
        count = (subset["P"] < threshold).sum()
        add(
            f"P-value < {threshold:g}: "
            f"{count:,} ({count / len(subset) * 100:.4f}%)"
        )


# ============================================================
# 4. HETEROZYGOSITY AND F
# ============================================================

het["OBS_HET"] = (
    (het["N(NM)"] - het["O(HOM)"]) / het["N(NM)"]
)

het = het.replace([np.inf, -np.inf], np.nan)
het = het.dropna(subset=["OBS_HET", "F"])

summarize("OBSERVED HETEROZYGOSITY", het["OBS_HET"])
summarize("INBREEDING COEFFICIENT (F)", het["F"])

f_mean = het["F"].mean()
f_std = het["F"].std()

lower = f_mean - 3 * f_std
upper = f_mean + 3 * f_std

outliers = het[(het["F"] < lower) | (het["F"] > upper)]

add("\nF OUTLIERS (MEAN +/- 3 SD)")
add(f"Mean F: {f_mean:.6f}")
add(f"SD F: {f_std:.6f}")
add(f"Lower bound: {lower:.6f}")
add(f"Upper bound: {upper:.6f}")
add(f"Number of outliers: {len(outliers):,}")

if not outliers.empty:
    add(outliers[["FID", "IID", "F", "OBS_HET"]].to_string(index=False))


# ============================================================
# PLOT 1: MAF DISTRIBUTION (LOG SCALE)
# ============================================================

maf_positive = maf[maf > 0]

plt.figure(figsize=(8, 5))
plt.hist(maf_positive, bins=np.logspace(
    np.log10(maf_positive.min()),
    np.log10(maf_positive.max()),
    60
), edgecolor="black")

plt.xscale("log")
plt.xlabel("Minor Allele Frequency (MAF)")
plt.ylabel("Number of Variants")
plt.title("MAF Distribution")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(PLOT_DIR / "01_maf_distribution.png", dpi=300)
plt.close()


# ============================================================
# PLOT 2: HWE P-VALUE DISTRIBUTION
# ============================================================

p_values = hwe["P"].dropna()
neg_log_p = -np.log10(p_values.clip(lower=1e-300))

plt.figure(figsize=(8, 5))
plt.hist(neg_log_p, bins=60, edgecolor="black")

plt.axvline(-np.log10(1e-6), linestyle="--",
            label="P = 1e-6")
plt.axvline(-np.log10(1e-10), linestyle=":",
            label="P = 1e-10")

plt.xlabel("-log10(HWE P-value)")
plt.ylabel("Number of Variants")
plt.title("HWE P-value Distribution")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(PLOT_DIR / "02_hwe_distribution.png", dpi=300)
plt.close()



# ============================================================
# PLOT 3: MISSINGNESS HISTOGRAMS
# ============================================================

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Sample Missingness
x = imiss["F_MISS"].dropna()

axes[0].hist(x, bins=50, edgecolor="black", alpha=0.75)
axes[0].axvline(x.mean(), color="red", linestyle="--",
                label=f"Mean = {x.mean():.6f}")
axes[0].axvline(x.median(), color="green", linestyle=":",
                label=f"Median = {x.median():.6f}")

axes[0].set_title("Sample Missingness")
axes[0].set_xlabel("Missingness Rate")
axes[0].set_ylabel("Number of Individuals")
axes[0].legend()
axes[0].grid(axis="y", alpha=0.3)


# Variant Missingness (zoomed)
x = lmiss["F_MISS"].dropna()

axes[1].hist(
    x,
    bins=100,
    range=(0, 0.05),
    edgecolor="black",
    alpha=0.75
)

for threshold, color in [
    (0.01, "red"),
    (0.02, "orange"),
    (0.05, "green")
]:
    count = (x > threshold).sum()
    percentage = count / len(x) * 100

    axes[1].axvline(
        threshold,
        color=color,
        linestyle="--",
        linewidth=2,
        label=f">{threshold:.0%}: {count:,} ({percentage:.4f}%)"
    )

axes[1].set_title("Variant Missingness (Zoomed: 0–5%)")
axes[1].set_xlabel("Missingness Rate")
axes[1].set_ylabel("Number of Variants")
axes[1].legend(fontsize=8)
axes[1].grid(axis="y", alpha=0.3)

plt.tight_layout()
plt.savefig(
    PLOT_DIR / "03_missingness_histograms.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()



# ============================================================
# PLOT 4: F VS OBSERVED HETEROZYGOSITY
# ============================================================

plt.figure(figsize=(8, 5))

plt.scatter(
    het["OBS_HET"],
    het["F"],
    alpha=0.5,
    s=15,
    label="Individuals"
)

plt.axhline(upper, linestyle="--", label="+3 SD")
plt.axhline(lower, linestyle="--", label="-3 SD")
plt.axhline(f_mean, linestyle=":", label="Mean F")

if not outliers.empty:
    plt.scatter(
        outliers["OBS_HET"],
        outliers["F"],
        color="red",
        s=25,
        label="F outliers"
    )

    for _, row in outliers.iterrows():
        plt.annotate(
            row["IID"],
            (row["OBS_HET"], row["F"]),
            fontsize=7
        )

plt.xlabel("Observed Heterozygosity (Ho)")
plt.ylabel("Inbreeding Coefficient (F)")
plt.title("F vs Observed Heterozygosity")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(PLOT_DIR / "04_f_vs_heterozygosity.png", dpi=300)
plt.close()


# ============================================================
# SAVE REPORT
# ============================================================

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text("\n".join(report) + "\n")

print("\n".join(report))
print(f"\nReport saved to: {OUTPUT}")
print(f"Plots saved to: {PLOT_DIR}")
