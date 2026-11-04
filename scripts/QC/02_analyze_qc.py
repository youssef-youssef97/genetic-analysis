import pandas as pd
from pathlib import Path

BASE = Path("results/QC/chr14_initial")
OUTPUT = Path("results/QC/chr14_qc_summary.txt")

# Load initial QC reports
imiss = pd.read_csv(f"{BASE}.imiss", sep=r"\s+")
lmiss = pd.read_csv(f"{BASE}.lmiss", sep=r"\s+")
frq = pd.read_csv(f"{BASE}.frq", sep=r"\s+")
hwe = pd.read_csv(f"{BASE}.hwe", sep=r"\s+")

report = []

def add(text=""):
    report.append(str(text))

def summarize_missingness(name, df):
    x = df["F_MISS"]

    add(f"\n{'=' * 60}")
    add(name)
    add(f"Number of records: {len(x):,}")
    add(x.describe(
        percentiles=[0.50, 0.90, 0.95, 0.99, 0.999]
    ).to_string())

    for threshold in [0.01, 0.02, 0.05]:
        count = (x > threshold).sum()
        add(
            f"Records with missingness > {threshold:.0%}: "
            f"{count:,} ({count / len(x) * 100:.4f}%)"
        )

# Sample missingness
summarize_missingness("SAMPLE MISSINGNESS", imiss)

# Variant missingness
summarize_missingness("VARIANT MISSINGNESS", lmiss)

# Minor allele frequency
maf = frq["MAF"]

add(f"\n{'=' * 60}")
add("MINOR ALLELE FREQUENCY")
add(f"Number of variants: {len(maf):,}")
add(maf.describe(
    percentiles=[0.01, 0.05, 0.25, 0.50, 0.75, 0.95]
).to_string())

for threshold in [0.01, 0.05]:
    count = (maf < threshold).sum()
    add(
        f"Variants with MAF < {threshold:.0%}: "
        f"{count:,} ({count / len(maf) * 100:.4f}%)"
    )

# Hardy-Weinberg equilibrium
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

# Save report
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text("\n".join(report) + "\n")

print("\n".join(report))
print(f"\nReport saved to: {OUTPUT}")