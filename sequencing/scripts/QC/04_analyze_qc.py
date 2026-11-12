
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from tqdm import tqdm

# ============================================================
# PATHS
# ============================================================

BASE = Path("results/QC/genomewide_initial")
OUTPUT = Path("results/QC/genomewide_qc_summary.txt")
PLOT_DIR = Path("results/QC/qc_plots")

PLOT_DIR.mkdir(parents=True, exist_ok=True)

# Overall progress across the main analysis stages
progress = tqdm(
    total=7,
    desc="Genome-wide QC",
    unit="stage",
    dynamic_ncols=True
)


# ============================================================
# LOAD PLINK 2.0 OUTPUTS
# ============================================================

smiss = pd.read_csv(f"{BASE}.smiss", sep=r"\s+")
smiss.rename(columns={"#FID": "FID"}, inplace=True)

vmiss = pd.read_csv(f"{BASE}.vmiss", sep=r"\s+")
vmiss.rename(columns={"#CHROM": "CHROM"}, inplace=True)

afreq = pd.read_csv(f"{BASE}.afreq", sep=r"\s+")
afreq.rename(columns={"#CHROM": "CHROM"}, inplace=True)

hardy = pd.read_csv(f"{BASE}.hardy", sep=r"\s+")
hardy.rename(columns={"#CHROM": "CHROM"}, inplace=True)

het = pd.read_csv(f"{BASE}.het", sep=r"\s+")
het.rename(columns={"#FID": "FID"}, inplace=True)

progress.update(1)
progress.set_postfix_str("Input files loaded")

report = []


def add(text=""):
    report.append(str(text))


# ============================================================
# VALIDATE COLUMN NAMES
# ============================================================

def check_columns(df, columns, filename):
    missing = [col for col in columns if col not in df.columns]

    if missing:
        raise ValueError(
            f"Missing columns in {filename}: {missing}\n"
            f"Available columns: {list(df.columns)}"
        )


check_columns(
    smiss,
    ["FID", "IID", "F_MISS"],
    "smiss"
)

check_columns(
    vmiss,
    ["CHROM", "ID", "F_MISS"],
    "vmiss"
)

check_columns(
    afreq,
    ["CHROM", "ID", "ALT_FREQS"],
    "afreq"
)

check_columns(
    hardy,
    ["CHROM", "ID", "P"],
    "hardy"
)

check_columns(
    het,
    ["FID", "IID", "O(HOM)", "OBS_CT", "F"],
    "het"
)
progress.update(1)
progress.set_postfix_str("Column validation complete")


# ============================================================
# GENERAL SUMMARY FUNCTION
# ============================================================

def summarize(name, x, thresholds=None, direction="greater"):

    x = pd.to_numeric(x, errors="coerce").dropna()

    add(f"\n{'=' * 60}")
    add(name)
    add(f"Number of records: {len(x):,}")

    if x.empty:
        add("No valid records available.")
        return

    add(
        x.describe(
            percentiles=[
                0.01, 0.05, 0.25, 0.50,
                0.75, 0.90, 0.95, 0.99
            ]
        ).to_string()
    )

    if thresholds:

        for t in thresholds:

            if direction == "less":
                count = (x < t).sum()
                operator = "<"
            else:
                count = (x > t).sum()
                operator = ">"

            percentage = count / len(x) * 100

            add(
                f"Records {operator} {t:g}: "
                f"{count:,} ({percentage:.4f}%)"
            )


# ============================================================
# 1. SAMPLE AND VARIANT MISSINGNESS
# ============================================================

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
progress.update(1)
progress.set_postfix_str("Missingness summary complete")


# ============================================================
# 2. MINOR ALLELE FREQUENCY (MAF)
# ============================================================

def calculate_maf(alt_freqs):

    """
    Calculate the frequency of the least frequent allele.

    ALT_FREQS contains comma-separated frequencies of
    all alternate alleles.

    REF frequency = 1 - sum(ALT frequencies).
    """

    if pd.isna(alt_freqs):
        return np.nan

    try:

        alt_values = [
            float(value)
            for value in str(alt_freqs).split(",")
        ]

        ref_freq = 1.0 - sum(alt_values)

        all_frequencies = [ref_freq] + alt_values

        if any(
            not np.isfinite(value)
            or value < -1e-6
            or value > 1 + 1e-6
            for value in all_frequencies
        ):
            return np.nan

        return min(all_frequencies)

    except (ValueError, TypeError):
        return np.nan


tqdm.pandas(desc="Calculating MAF", dynamic_ncols=True)
maf = afreq["ALT_FREQS"].progress_apply(calculate_maf)
maf = maf.dropna()

summarize(
    "MINOR ALLELE FREQUENCY",
    maf,
    [0.01, 0.05],
    direction="less"
)
progress.update(1)
progress.set_postfix_str("MAF analysis complete")


# ============================================================
# 3. HARDY-WEINBERG EQUILIBRIUM
# ============================================================

add(f"\n{'=' * 60}")
add("HARDY-WEINBERG EQUILIBRIUM")

add(f"Number of records: {len(hardy):,}")

hwe_p = pd.to_numeric(
    hardy["P"],
    errors="coerce"
).dropna()

add("\nP-value distribution:")
add(hwe_p.describe().to_string())

for threshold in [1e-6, 1e-10]:

    count = (hwe_p < threshold).sum()

    percentage = (
        count / len(hwe_p) * 100
        if len(hwe_p) else 0
    )

    add(
        f"P-value < {threshold:g}: "
        f"{count:,} ({percentage:.4f}%)"
    )


# HWE summary per chromosome

add("\nHWE P-value summary per chromosome:")

hwe_summary = hardy.groupby("CHROM")["P"].agg(
    ["count", "median", "min"]
)

add(hwe_summary.to_string())
progress.update(1)
progress.set_postfix_str("HWE analysis complete")


# ============================================================
# 4. HETEROZYGOSITY AND INBREEDING COEFFICIENT
# ============================================================

het["OBS_HET"] = (
    pd.to_numeric(het["OBS_CT"], errors="coerce")
    -
    pd.to_numeric(het["O(HOM)"], errors="coerce")
) / pd.to_numeric(het["OBS_CT"], errors="coerce")

het["F"] = pd.to_numeric(
    het["F"],
    errors="coerce"
)

het = het.replace(
    [np.inf, -np.inf],
    np.nan
)

het = het.dropna(
    subset=["OBS_HET", "F"]
)

summarize(
    "OBSERVED HETEROZYGOSITY",
    het["OBS_HET"]
)

summarize(
    "INBREEDING COEFFICIENT (F)",
    het["F"]
)


# F outlier detection: Mean +/- 3 SD

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
        outliers[
            ["FID", "IID", "F", "OBS_HET"]
        ].to_string(index=False)
    )

progress.update(1)
progress.set_postfix_str("Heterozygosity analysis complete")


# ============================================================
# PLOT 1: MAF DISTRIBUTION (LOG SCALE)
# ============================================================

maf_positive = maf[maf > 0]

if not maf_positive.empty:

    plt.figure(figsize=(8, 5))

    plt.hist(
        maf_positive,
        bins=np.logspace(
            np.log10(maf_positive.min()),
            np.log10(maf_positive.max()),
            60
        ),
        edgecolor="black"
    )

    plt.xscale("log")
    plt.xlabel("Minor Allele Frequency (MAF)")
    plt.ylabel("Number of Variants")
    plt.title("MAF Distribution")
    plt.grid(alpha=0.3)
    plt.tight_layout()

    plt.savefig(
        PLOT_DIR / "01_maf_distribution.png",
        dpi=300
    )

    plt.close()


# ============================================================
# PLOT 2: HWE P-VALUE DISTRIBUTION
# ============================================================

p_values = hwe_p

if not p_values.empty:

    neg_log_p = -np.log10(
        p_values.clip(lower=1e-300)
    )

    plt.figure(figsize=(8, 5))

    plt.hist(
        neg_log_p,
        bins=60,
        edgecolor="black"
    )

    plt.axvline(
        -np.log10(1e-6),
        linestyle="--",
        label="P = 1e-6"
    )

    plt.axvline(
        -np.log10(1e-10),
        linestyle=":",
        label="P = 1e-10"
    )

    plt.xlabel("-log10(HWE P-value)")
    plt.ylabel("Number of Variants")
    plt.title("HWE P-value Distribution")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()

    plt.savefig(
        PLOT_DIR / "02_hwe_distribution.png",
        dpi=300
    )

    plt.close()


# ============================================================
# PLOT 3: MISSINGNESS HISTOGRAMS
# ============================================================

fig, axes = plt.subplots(
    1, 2,
    figsize=(14, 5)
)


# Sample Missingness

x = smiss["F_MISS"].dropna()

axes[0].hist(
    x,
    bins=50,
    edgecolor="black",
    alpha=0.75
)

if not x.empty:

    axes[0].axvline(
        x.mean(),
        linestyle="--",
        label=f"Mean = {x.mean():.6f}"
    )

    axes[0].axvline(
        x.median(),
        linestyle=":",
        label=f"Median = {x.median():.6f}"
    )

axes[0].set_title("Sample Missingness")
axes[0].set_xlabel("Missingness Rate")
axes[0].set_ylabel("Number of Individuals")
axes[0].legend()
axes[0].grid(axis="y", alpha=0.3)


# Variant Missingness (Zoomed)

x = vmiss["F_MISS"].dropna()

axes[1].hist(
    x,
    bins=100,
    range=(0, 0.05),
    edgecolor="black",
    alpha=0.75
)

for threshold in [0.01, 0.02, 0.05]:

    count = (x > threshold).sum()

    percentage = (
        count / len(x) * 100
        if len(x) else 0
    )

    axes[1].axvline(
        threshold,
        linestyle="--",
        linewidth=2,
        label=(
            f">{threshold:.0%}: "
            f"{count:,} ({percentage:.4f}%)"
        )
    )

axes[1].set_title(
    "Variant Missingness (Zoomed: 0–5%)"
)

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

plt.axhline(
    upper,
    linestyle="--",
    label="+3 SD"
)

plt.axhline(
    lower,
    linestyle="--",
    label="-3 SD"
)

plt.axhline(
    f_mean,
    linestyle=":",
    label="Mean F"
)

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

plt.savefig(
    PLOT_DIR / "04_f_vs_heterozygosity.png",
    dpi=300
)

plt.close()


# ============================================================
# SAVE REPORT
# ============================================================

OUTPUT.parent.mkdir(
    parents=True,
    exist_ok=True
)

OUTPUT.write_text(
    "\n".join(report) + "\n"
)

print("\n".join(report))

print(f"\nReport saved to: {OUTPUT}")
print(f"Plots saved to: {PLOT_DIR}")

progress.update(1)
progress.set_postfix_str("Report and plots saved")
progress.close()
