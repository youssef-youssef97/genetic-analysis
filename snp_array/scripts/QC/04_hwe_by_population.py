
import pandas as pd
import subprocess
from pathlib import Path

# ==================================================
# PATHS
# ==================================================

PROJECT = Path("snp_array")

PLINK2 = "/home/youssef/Project/plink2_install/plink2"

GENOTYPE = PROJECT / "data/processed/matched/1kg_omni_matched"
METADATA = PROJECT / "data/raw/plink/metadata.psam"

RESULTS = PROJECT / "results/QC/HWE"

SUPER_DIR = RESULTS / "superpopulation"
POP_DIR = RESULTS / "population"

SUPER_DIR.mkdir(parents=True, exist_ok=True)
POP_DIR.mkdir(parents=True, exist_ok=True)

# ==================================================
# LOAD SAMPLES AND METADATA
# ==================================================

print("Loading sample information...")

fam = pd.read_csv(
    f"{GENOTYPE}.fam",
    sep=r"\s+",
    header=None,
    names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"]
)

metadata = pd.read_csv(
    METADATA,
    sep=r"\s+"
).rename(columns={"#IID": "IID"})

samples = fam[["FID", "IID"]].merge(
    metadata[["IID", "SuperPop", "Population"]],
    on="IID",
    how="inner",
    validate="one_to_one"
)

print("Matched samples:", len(samples))

# ==================================================
# FUNCTION: RUN HWE FOR EACH GROUP
# ==================================================

summary = []


def run_hwe(group_column, output_dir):

    groups = sorted(samples[group_column].dropna().unique())

    print(f"\nRunning HWE by {group_column}...")

    for group in groups:

        group_samples = samples[
            samples[group_column] == group
        ].copy()

        keep_file = output_dir / f"{group}.keep"
        output_prefix = output_dir / group

        # PLINK keep file: FID IID
        group_samples[["FID", "IID"]].to_csv(
            keep_file,
            sep="\t",
            index=False,
            header=False
        )

        print(
            f"{group}: {len(group_samples)} samples"
        )

        # Run PLINK 2 HWE test
        command = [
            PLINK2,
            "--bfile", str(GENOTYPE),
            "--keep", str(keep_file),
            "--hardy",
            "--out", str(output_prefix)
        ]

        subprocess.run(command, check=True)

        # Read HWE results
        hardy = pd.read_csv(
            f"{output_prefix}.hardy",
            sep=r"\s+"
        )

        p_values = pd.to_numeric(
            hardy["P"],
            errors="coerce"
        ).dropna()

        
        count_1e6 = (p_values < 1e-6).sum()
        count_1e10 = (p_values < 1e-10).sum()
        total = len(p_values)

        summary.append({
            "Level": group_column,
            "Group": group,
            "Samples": len(group_samples),
            "Variants_tested": total,
            "P_below_1e-6": count_1e6,
            "P_below_1e-6_pct": (
                count_1e6 / total * 100 if total else 0
            ),
            "P_below_1e-10": count_1e10,
            "P_below_1e-10_pct": (
                count_1e10 / total * 100 if total else 0
            )
        })


        print(
            f"  Variants tested: {len(p_values):,}"
        )

        print(
            f"  P < 1e-6: {(p_values < 1e-6).sum():,}"
        )


# ==================================================
# 1. HWE BY SUPERPOPULATION
# ==================================================

run_hwe("SuperPop", SUPER_DIR)


# ==================================================
# 2. HWE BY POPULATION
# ==================================================

run_hwe("Population", POP_DIR)


# ==================================================
# SAVE SUMMARY
# ==================================================

summary_df = pd.DataFrame(summary)

summary_file = RESULTS / "HWE_summary.tsv"

summary_df.to_csv(
    summary_file,
    sep="\t",
    index=False
)

print("\n======================================")
print("HWE ANALYSIS COMPLETED")
print("======================================")

print(summary_df.to_string(index=False))

print("\nSummary saved:", summary_file)
print("Results directory:", RESULTS)
