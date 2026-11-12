
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# PATHS
# ============================================================

BASE = Path("snp_array/results/QC/1kg_omni_auto")

OUTPUT = Path("snp_array/results/QC/qc_summary.txt")
PLOT_DIR = Path("snp_array/results/QC/qc_plots")

PLOT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# LOAD PLINK 2.0 OUTPUTS
# ============================================================

print("Loading PLINK QC results...")

smiss = pd.read_csv(f"{BASE}.smiss", sep=r"\s+")
vmiss = pd.read_csv(f"{BASE}.vmiss", sep=r"\s+")
afreq = pd.read_csv(f"{BASE}.afreq", sep=r"\s+")
hardy = pd.read_csv(f"{BASE}.hardy", sep=r"\s+")
het = pd.read_csv(f"{BASE}.het", sep=r"\s+")

# Normalize column names
for df in [smiss, het]:
    df.rename(columns={"#FID": "FID"}, inplace=True)

for df in [vmiss, afreq, hardy]:
    df.rename(columns={"#CHROM": "CHROM"}, inplace=True)

report = []

def add(text=""):
    report.append(str(text))


# ============================================================
# GENERAL SUMMARY
# ============================================================

def summarize(name, values, thresholds=None, direction="greater"):

    x = pd.to_numeric(values, errors="coerce").dropna()

    add("\n" + "=" * 60)
    add(name)
    add(f"Number of records: {len(x):,}")

    if x.empty:
        add("No valid records.")
        return

    add(x.describe(
        percentiles=[0.01, 0.05, 0.25, 0.50,
                     0.75, 0.90, 0.95, 0.99]
    ).to_string())

    if thresholds:
        for t in thresholds:

            count = (x < t).sum() if direction == "less" else (x > t).sum()
            operator = "<" if direction == "less" else ">"

            add(
                f"Records {operator} {t:g}: "
                f"{count:,} ({count / len(x) * 100:.4f}%)"
            )


# ============================================================
# 1. MISSINGNESS
# ============================================================

print("Analyzing missingness...")

summarize(
    "SAMPLE MISSINGNESS",
    smiss["F_MISS"],
    [0.01, 0.02, 0.05]
)

summarize(
    "VARIANT MISSINGNESS",
    vmiss["F_MISS"],
    [0.01, 0.02, 0.05]
)


# ============================================================
# 2. MINOR ALLELE FREQUENCY
# ============================================================

print("Analyzing allele frequencies...")

def calculate_maf(value):

    if pd.isna(value):
        return np.nan

    try:
        alt = [float(v) for v in str(value).split(",")]
        ref = 1.0 - sum(alt)

        frequencies = [ref] + alt

        if any(v < -1e-6 or v > 1 + 1e-6 for v in frequencies):
            return np.nan

        return min(frequencies)

    except (ValueError, TypeError):
        return np.nan


maf = afreq["ALT_FREQS"].apply(calculate_maf).dropna()

summarize(
    "MINOR ALLELE FREQUENCY",
    maf,
    [0.01, 0.05],
    direction="less"
)


# ============================================================
# 3. HARDY-WEINBERG EQUILIBRIUM
# ============================================================

print("Analyzing Hardy-Weinberg equilibrium...")

hwe_p = pd.to_numeric(hardy["P"], errors="coerce").dropna()

add("\n" + "=" * 60)
add("HARDY-WEINBERG EQUILIBRIUM")
add(f"Number of variants: {len(hwe_p):,}")
add(hwe_p.describe().to_string())

for threshold in [1e-6, 1e-10]:

    count = (hwe_p < threshold).sum()

    add(
        f"P-value < {threshold:g}: "
        f"{count:,} ({count / len(hwe_p) * 100:.4f}%)"
        if len(hwe_p) else f"P-value < {threshold:g}: 0"
    )


# ============================================================
# 4. HETEROZYGOSITY AND INBREEDING COEFFICIENT
# ============================================================

print("Analyzing heterozygosity...")

het["OBS_CT"] = pd.to_numeric(het["OBS_CT"], errors="coerce")
het["O(HOM)"] = pd.to_numeric(het["O(HOM)"], errors="coerce")
het["F"] = pd.to_numeric(het["F"], errors="coerce")

het["OBS_HET"] = (
    het["OBS_CT"] - het["O(HOM)"]
) / het["OBS_CT"]

het = het.replace([np.inf, -np.inf], np.nan)
het = het.dropna(subset=["OBS_HET", "F"])

summarize(
    "OBSERVED HETEROZYGOSITY",
    het["OBS_HET"]
)

summarize(
    "INBREEDING COEFFICIENT (F)",
    het["F"]
)

# F outliers: Mean +/- 3 SD

f_mean = het["F"].mean()
f_std = het["F"].std()

lower = f_mean - 3 * f_std
upper = f_mean + 3 * f_std

outliers = het[
    (het["F"] < lower) |
    (het["F"] > upper)
]

add("\nF OUTLIERS (MEAN +/- 3 SD)")
add(f"Mean F: {f_mean:.6f}")
add(f"SD F: {f_std:.6f}")
add(f"Lower bound: {lower:.6f}")
add(f"Upper bound: {upper:.6f}")
add(f"Number of outliers: {len(outliers):,}")

if not outliers.empty:
    add(
        outliers[["FID", "IID", "F", "OBS_HET"]]
        .to_string(index=False)
    )


# ============================================================
# PLOT 1: MAF DISTRIBUTION
# ============================================================

maf_positive = maf[maf > 0]

if not maf_positive.empty:

    plt.figure(figsize=(8, 5))

    plt.hist(
        maf_positive,
        bins=60,
        edgecolor="black"
    )

    plt.xlabel("Minor Allele Frequency")
    plt.ylabel("Number of Variants")
    plt.title("MAF Distribution")
    plt.tight_layout()

    plt.savefig(PLOT_DIR / "01_maf_distribution.png", dpi=300)
    plt.close()


# ============================================================
# PLOT 2: HWE DISTRIBUTION
# ============================================================

if not hwe_p.empty:

    plt.figure(figsize=(8, 5))

    plt.hist(
        -np.log10(hwe_p.clip(lower=1e-300)),
        bins=60,
        edgecolor="black"
    )

    plt.axvline(-np.log10(1e-6), linestyle="--", label="P = 1e-6")
    plt.axvline(-np.log10(1e-10), linestyle=":", label="P = 1e-10")

    plt.xlabel("-log10(HWE P-value)")
    plt.ylabel("Number of Variants")
    plt.title("HWE P-value Distribution")
    plt.legend()
    plt.tight_layout()

    plt.savefig(PLOT_DIR / "02_hwe_distribution.png", dpi=300)
    plt.close()


# ============================================================
# PLOT 3: MISSINGNESS
# ============================================================

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Sample missingness

axes[0].hist(smiss["F_MISS"].dropna(), bins=50, edgecolor="black")
axes[0].set_title("Sample Missingness")
axes[0].set_xlabel("Missingness Rate")
axes[0].set_ylabel("Number of Individuals")

# Variant missingness

axes[1].hist(
    vmiss["F_MISS"].dropna(),
    bins=100,
    range=(0, 0.05),
    edgecolor="black"
)

axes[1].set_title("Variant Missingness (0-5%)")
axes[1].set_xlabel("Missingness Rate")
axes[1].set_ylabel("Number of Variants")

plt.tight_layout()
plt.savefig(PLOT_DIR / "03_missingness.png", dpi=300)
plt.close()


# ============================================================
# PLOT 4: F VS OBSERVED HETEROZYGOSITY
# ============================================================

plt.figure(figsize=(8, 5))

plt.scatter(
    het["OBS_HET"],
    het["F"],
    alpha=0.5,
    s=15
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
        label="Outliers"
    )

plt.xlabel("Observed Heterozygosity")
plt.ylabel("Inbreeding Coefficient (F)")
plt.title("F vs Observed Heterozygosity")
plt.legend()
plt.tight_layout()

plt.savefig(PLOT_DIR / "04_heterozygosity.png", dpi=300)
plt.close()


# ============================================================
# SAVE REPORT
# ============================================================

OUTPUT.write_text("\n".join(report) + "\n")

print("\n" + "\n".join(report))

print("\nQC analysis completed successfully!")
print(f"Report: {OUTPUT}")
print(f"Plots: {PLOT_DIR}")
