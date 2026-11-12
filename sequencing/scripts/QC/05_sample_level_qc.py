#!/usr/bin/env python3
"""Robust sample-level QC pipeline for 1000 Genomes Phase 3 with PLINK 2.

Stages: LD pruning, autosomal heterozygosity, KING relatedness,
chromosome-X sex concordance, chromosome-Y missingness, summary, plots.

Heterozygosity outliers are screening flags, not proof of contamination.
Population-specific MAD is calculated when the 1000G panel is available.
The Y fileset is male-only, so its missingness is descriptive QC, not an
independent sex-concordance test.
"""
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from tqdm import tqdm


# --------------------------- CONFIGURATION ---------------------------
PROJECT = Path("/home/youssef/Project/genetic-analysis")
PLINK2 = Path("/home/youssef/Project/plink2_install/plink2")

AUTO = PROJECT / "data/raw/chromosomes/plink/autosomes/autosomes"
CHR_X = PROJECT / "data/raw/chromosomes/plink/chrX"
CHR_Y = PROJECT / "data/raw/chromosomes/plink/chrY"
PANEL = PROJECT / (
    "data/raw/chromosomes/vcf/"
    "integrated_call_samples_v3.20130502.ALL.panel"
)

OUTPUT_DIR = PROJECT / "results/QC/sample_level"
PLOT_DIR = OUTPUT_DIR / "plots"
LOG_DIR = OUTPUT_DIR / "logs"
REPORT_FILE = OUTPUT_DIR / "sample_level_qc_report.txt"

PRUNE_WINDOW = 200
PRUNE_STEP = 50
PRUNE_R2 = 0.1
MIN_MAF = 0.01

KING_DUPLICATE = 0.354
KING_FIRST_DEGREE = 0.177
KING_SECOND_DEGREE = 0.0884
KING_THIRD_DEGREE = 0.0442

HET_MAD_THRESHOLD = 3.0
MAX_FEMALE_XF = 0.2
MIN_MALE_XF = 0.8
Y_MISSINGNESS_THRESHOLD = 0.05

for directory in (OUTPUT_DIR, PLOT_DIR, LOG_DIR):
    directory.mkdir(parents=True, exist_ok=True)


# --------------------------- HELPERS ---------------------------
report = []


def add(message=""):
    report.append(str(message))


def section(title):
    add()
    add("=" * 78)
    add(title)
    add("=" * 78)


def save_report():
    REPORT_FILE.write_text("\n".join(report) + "\n", encoding="utf-8")


def require_file(path):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Required file not found: {path}")
    return path


def require_pfile(prefix):
    missing = [
        str(prefix) + ext for ext in (".pgen", ".pvar", ".psam")
        if not Path(str(prefix) + ext).exists()
    ]
    if missing:
        raise FileNotFoundError(
            "Incomplete PLINK 2 fileset. Missing:\n  " + "\n  ".join(missing)
        )


def run_plink(command, stage_name):
    log_file = LOG_DIR / f"{stage_name}.log"
    tqdm.write(f"Running: {stage_name}")
    with log_file.open("w", encoding="utf-8") as log:
        result = subprocess.run(
            [str(x) for x in command],
            cwd=str(PROJECT),
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
    if result.returncode != 0:
        lines = log_file.read_text(
            encoding="utf-8", errors="replace"
        ).splitlines()
        raise RuntimeError(
            f"PLINK 2 failed in {stage_name} (exit {result.returncode}).\n"
            f"Log: {log_file}\n\n" + "\n".join(lines[-35:])
        )
    tqdm.write(f"Completed: {stage_name}")


def read_plink_table(path):
    """Read PLINK whitespace table and normalize #FID/#CHROM headers."""
    frame = pd.read_csv(path, sep=r"\s+", dtype=str)
    frame.columns = [
        str(column).lstrip("#").strip().upper()
        for column in frame.columns
    ]
    return frame


def require_columns(frame, required, filename):
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(
            f"Unexpected columns in {filename}.\n"
            f"Missing: {missing}\nAvailable: {list(frame.columns)}"
        )


def to_numeric(frame, column):
    frame[column] = pd.to_numeric(frame[column], errors="coerce")


def robust_mad(values, threshold=3.0):
    """Median/MAD outlier flags; zero MAD yields no automatic flags."""
    values = pd.to_numeric(values, errors="coerce")
    valid = values.dropna()
    if valid.empty:
        flags = pd.Series(False, index=values.index)
        return np.nan, np.nan, np.nan, np.nan, np.nan, flags

    median = float(valid.median())
    mad = float(np.median(np.abs(valid.to_numpy() - median)))
    robust_sd = 1.4826 * mad

    if not np.isfinite(robust_sd) or robust_sd <= 0:
        lower = upper = median
        flags = pd.Series(False, index=values.index)
    else:
        lower = median - threshold * robust_sd
        upper = median + threshold * robust_sd
        flags = ((values < lower) | (values > upper)).fillna(False)

    return median, mad, robust_sd, lower, upper, flags


def classify_relationship(value):
    if value >= KING_DUPLICATE:
        return "Duplicate / MZ twin"
    if value >= KING_FIRST_DEGREE:
        return "First-degree"
    if value >= KING_SECOND_DEGREE:
        return "Second-degree"
    if value >= KING_THIRD_DEGREE:
        return "Third-degree"
    return "Below third-degree threshold"


# --------------------------- INPUT VALIDATION ---------------------------
require_file(PLINK2)
require_pfile(AUTO)
require_pfile(CHR_X)
require_pfile(CHR_Y)

stages = [
    "LD pruning", "Autosomal heterozygosity", "KING relatedness",
    "Sex consistency (chrX)", "Y chromosome sample QC",
    "Statistical analysis", "Report and plots"
]
progress = tqdm(total=len(stages), desc="Sample-level QC", unit="stage")

try:
    # --------------------------- STAGE 1 ---------------------------
    section("STAGE 1: LD PRUNING")
    prune_prefix = OUTPUT_DIR / "autosomal_pruned"
    run_plink([
        PLINK2, "--pfile", AUTO,
        "--snps-only", "just-acgt",
        "--max-alleles", "2",
        "--maf", str(MIN_MAF),
        "--indep-pairwise", str(PRUNE_WINDOW),
        str(PRUNE_STEP), str(PRUNE_R2),
        "--out", prune_prefix,
    ], "01_ld_pruning")

    prune_file = require_file(Path(str(prune_prefix) + ".prune.in"))
    with prune_file.open(encoding="utf-8") as handle:
        pruned_count = sum(1 for line in handle if line.strip())
    if pruned_count == 0:
        raise RuntimeError("The .prune.in file is empty; check input/filter settings.")
    add(f"LD-pruned SNP count: {pruned_count:,}")
    progress.update(1)

    # --------------------------- STAGE 2 ---------------------------
    section("STAGE 2: AUTOSOMAL HETEROZYGOSITY")
    het_prefix = OUTPUT_DIR / "autosomal_heterozygosity"
    run_plink([
        PLINK2, "--pfile", AUTO,
        "--extract", prune_file,
        "--het", "--out", het_prefix,
    ], "02_autosomal_heterozygosity")

    het_file = require_file(Path(str(het_prefix) + ".het"))
    het = read_plink_table(het_file)
    require_columns(het, ["FID", "IID", "O(HOM)", "OBS_CT", "F"], het_file.name)
    for col in ("O(HOM)", "OBS_CT", "F"):
        to_numeric(het, col)

    het["OBS_HET"] = np.where(
        het["OBS_CT"] > 0,
        (het["OBS_CT"] - het["O(HOM)"]) / het["OBS_CT"],
        np.nan,
    )
    het.replace([np.inf, -np.inf], np.nan, inplace=True)
    valid_het = het.dropna(subset=["OBS_HET", "F"]).copy()
    if valid_het.empty:
        raise RuntimeError("No valid sample heterozygosity values in .het.")

    (het_median, het_mad, het_robust_sd, het_lower,
     het_upper, global_flags) = robust_mad(valid_het["OBS_HET"],
                                            HET_MAD_THRESHOLD)
    valid_het["GLOBAL_HET_OUTLIER"] = global_flags.to_numpy()
    global_outliers = valid_het[valid_het["GLOBAL_HET_OUTLIER"]].copy()

    # Population-specific analysis, if the panel file is present and valid.
    population_available = False
    if PANEL.exists():
        panel = pd.read_csv(PANEL, sep=r"\s+", dtype=str)
        panel.columns = [str(c).strip().upper() for c in panel.columns]
        if {"SAMPLE", "POP"}.issubset(panel.columns):
            panel = panel[["SAMPLE", "POP"]].drop_duplicates("SAMPLE")
            panel.rename(columns={"SAMPLE": "IID", "POP": "POPULATION"},
                         inplace=True)
            valid_het = valid_het.merge(panel, on="IID", how="left")
            population_available = True
            valid_het["POP_HET_OUTLIER"] = False
            valid_het["POP_HET_MEDIAN"] = np.nan
            valid_het["POP_HET_MAD"] = np.nan

            for pop, indices in valid_het.groupby(
                "POPULATION", dropna=True
            ).groups.items():
                _, pop_mad, _, _, _, pop_flags = robust_mad(
                    valid_het.loc[indices, "OBS_HET"], HET_MAD_THRESHOLD
                )
                pop_values = valid_het.loc[indices, "OBS_HET"]
                valid_het.loc[indices, "POP_HET_MEDIAN"] = pop_values.median()
                valid_het.loc[indices, "POP_HET_MAD"] = pop_mad
                valid_het.loc[indices, "POP_HET_OUTLIER"] = pop_flags.to_numpy()

            population_outliers = valid_het[
                valid_het["POP_HET_OUTLIER"]
            ].copy()
            add("1000 Genomes panel found: population-specific MAD was calculated.")
            add("Population-specific outlier counts:")
            if population_outliers.empty:
                add("No population-specific heterozygosity outliers.")
            else:
                add(population_outliers["POPULATION"].value_counts().to_string())
            add(f"Samples without panel population: "
                f"{valid_het['POPULATION'].isna().sum():,}")
            population_outliers.to_csv(
                OUTPUT_DIR / "heterozygosity_population_outliers.tsv",
                sep="\t", index=False
            )
        else:
            add("WARNING: panel file lacks SAMPLE/POP columns; population analysis skipped.")
    else:
        add("NOTE: panel file not found; population-specific analysis skipped.")

    add(f"Samples with valid heterozygosity: {len(valid_het):,}")
    add(f"Pooled median observed heterozygosity: {het_median:.8f}")
    add(f"Pooled MAD: {het_mad:.8f}")
    add(f"Pooled robust SD: {het_robust_sd:.8f}")
    add(f"Pooled lower threshold: {het_lower:.8f}")
    add(f"Pooled upper threshold: {het_upper:.8f}")
    add(f"Pooled heterozygosity outliers: {len(global_outliers):,}")

    global_outliers.to_csv(
        OUTPUT_DIR / "heterozygosity_global_outliers.tsv",
        sep="\t", index=False
    )
    valid_het.to_csv(
        OUTPUT_DIR / "heterozygosity_all_samples.tsv",
        sep="\t", index=False
    )
    progress.update(1)

    # --------------------------- STAGE 3 ---------------------------
    section("STAGE 3: DUPLICATION AND RELATEDNESS (KING)")
    king_prefix = OUTPUT_DIR / "king_relatedness"
    run_plink([
        PLINK2, "--pfile", AUTO,
        "--extract", prune_file,
        "--make-king-table",
        "--out", king_prefix,
    ], "03_king_relatedness")

    king_file = require_file(Path(str(king_prefix) + ".kin0"))
    king = read_plink_table(king_file)
    require_columns(king, ["KINSHIP"], king_file.name)
    to_numeric(king, "KINSHIP")
    king = king.dropna(subset=["KINSHIP"]).copy()
    king["RELATIONSHIP"] = king["KINSHIP"].apply(classify_relationship)

    related = king[king["KINSHIP"] >= KING_THIRD_DEGREE].copy()
    related.sort_values("KINSHIP", ascending=False, inplace=True)
    add(f"Total pairs in KING output: {len(king):,}")
    add(f"Pairs with kinship >= {KING_THIRD_DEGREE}: {len(related):,}")
    add("Relationship categories:")
    add(related["RELATIONSHIP"].value_counts().to_string()
        if not related.empty else "No pairs met the third-degree threshold.")

    related.to_csv(OUTPUT_DIR / "related_sample_pairs.tsv",
                   sep="\t", index=False)
    id_columns = [c for c in ("FID1", "IID1", "FID2", "IID2")
                  if c in related.columns]
    add("Top 20 strongest related pairs:")
    if related.empty or not id_columns:
        add("No reportable related pairs or pair ID columns.")
    else:
        add(related.head(20)[id_columns + [
            "KINSHIP", "RELATIONSHIP"
        ]].to_string(index=False))
    progress.update(1)

    # --------------------------- STAGE 4 ---------------------------
    section("STAGE 4: SEX CONSISTENCY (CHR X)")
    sex_prefix = OUTPUT_DIR / "sex_check"
    run_plink([
        PLINK2, "--pfile", CHR_X,
        "--check-sex",
        f"max-female-xf={MAX_FEMALE_XF}",
        f"min-male-xf={MIN_MALE_XF}",
        "--out", sex_prefix,
    ], "04_sex_check_chrX")

    sex_file = require_file(Path(str(sex_prefix) + ".sexcheck"))
    sex = read_plink_table(sex_file)
    require_columns(sex, ["FID", "IID", "PEDSEX", "SNPSEX", "STATUS", "F"],
                    sex_file.name)
    to_numeric(sex, "F")
    sex["STATUS"] = sex["STATUS"].astype(str).str.upper()
    sex_problems = sex[sex["STATUS"] != "OK"].copy()

    add(f"Samples checked: {len(sex):,}")
    add(f"Sex-check problems / ambiguous calls: {len(sex_problems):,}")
    add("Sex-check status:")
    add(sex["STATUS"].value_counts(dropna=False).to_string())
    add("Reported sex (PEDSEX):")
    add(sex["PEDSEX"].value_counts(dropna=False).to_string())
    add("Inferred sex (SNPSEX):")
    add(sex["SNPSEX"].value_counts(dropna=False).to_string())

    sex_problems.to_csv(OUTPUT_DIR / "sex_check_problems.tsv",
                        sep="\t", index=False)
    sex.to_csv(OUTPUT_DIR / "sex_check_all_samples.tsv",
               sep="\t", index=False)
    progress.update(1)

    # --------------------------- STAGE 5 ---------------------------
    section("STAGE 5: CHROMOSOME Y SAMPLE QC")
    y_prefix = OUTPUT_DIR / "chrY_sample_missingness"
    run_plink([
        PLINK2, "--pfile", CHR_Y,
        "--missing", "sample-only",
        "--out", y_prefix,
    ], "05_chrY_sample_qc")

    y_file = require_file(Path(str(y_prefix) + ".smiss"))
    y_smiss = read_plink_table(y_file)
    require_columns(y_smiss, ["IID", "F_MISS"], y_file.name)
    to_numeric(y_smiss, "F_MISS")
    y_smiss = y_smiss.dropna(subset=["F_MISS"]).copy()

    add(f"Samples in male-only chrY fileset: {len(y_smiss):,}")
    if not y_smiss.empty:
        add(f"Mean Y missingness: {y_smiss['F_MISS'].mean():.8f}")
        add(f"Median Y missingness: {y_smiss['F_MISS'].median():.8f}")
        add(f"Samples with Y missingness > {Y_MISSINGNESS_THRESHOLD:.0%}: "
            f"{(y_smiss['F_MISS'] > Y_MISSINGNESS_THRESHOLD).sum():,}")
    else:
        add("WARNING: no valid Y missingness values.")
    add("Y missingness is descriptive QC among the male-only Y fileset; "
        "it is not an independent sex-concordance test.")

    y_smiss.to_csv(OUTPUT_DIR / "chrY_sample_missingness.tsv",
                   sep="\t", index=False)
    progress.update(1)

    # --------------------------- STAGE 6 ---------------------------
    section("STAGE 6: SAMPLE-LEVEL QC SUMMARY")
    add("1. AUTOSOMAL HETEROZYGOSITY")
    add(f"Samples analyzed: {len(valid_het):,}")
    add(f"Pooled outliers: {len(global_outliers):,}")
    if population_available:
        add(f"Population-specific outliers: {len(population_outliers):,}")

    add("\n2. KING RELATEDNESS")
    add(f"Duplicate/MZ threshold or higher: "
        f"{(related['KINSHIP'] >= KING_DUPLICATE).sum():,}")
    add(f"First-degree range: "
        f"{((related['KINSHIP'] >= KING_FIRST_DEGREE) & (related['KINSHIP'] < KING_DUPLICATE)).sum():,}")
    add(f"Second-degree range: "
        f"{((related['KINSHIP'] >= KING_SECOND_DEGREE) & (related['KINSHIP'] < KING_FIRST_DEGREE)).sum():,}")
    add(f"Third-degree range: "
        f"{((related['KINSHIP'] >= KING_THIRD_DEGREE) & (related['KINSHIP'] < KING_SECOND_DEGREE)).sum():,}")

    add("\n3. SEX CONSISTENCY")
    add(f"Sex-check problems / ambiguous calls: {len(sex_problems):,}")
    add("\n4. CHROMOSOME Y")
    add(f"Male samples in Y fileset: {len(y_smiss):,}")
    add(f"Y missingness > {Y_MISSINGNESS_THRESHOLD:.0%}: "
        f"{(y_smiss['F_MISS'] > Y_MISSINGNESS_THRESHOLD).sum():,}")

    add("\nINTERPRETATION CAUTIONS")
    add("- Heterozygosity outliers are review flags, not proof of contamination.")
    add("- Population-specific comparisons are preferable in this multi-ancestry dataset.")
    add("- KING relationship categories are approximate and need population context.")
    add("- Do not remove samples automatically based only on these flags.")
    progress.update(1)

    # --------------------------- STAGE 7: PLOTS ---------------------------
    section("STAGE 7: PLOTS")

    plt.figure(figsize=(9, 5))
    plt.hist(valid_het["OBS_HET"].dropna(), bins=50,
             edgecolor="black", alpha=0.75)
    plt.axvline(het_median, linestyle="--",
                label=f"Pooled median = {het_median:.5f}")
    if np.isfinite(het_lower) and np.isfinite(het_upper):
        plt.axvline(het_lower, linestyle=":", label="Pooled lower threshold")
        plt.axvline(het_upper, linestyle=":", label="Pooled upper threshold")
    plt.xlabel("Observed autosomal heterozygosity")
    plt.ylabel("Number of individuals")
    plt.title("Autosomal heterozygosity distribution")
    plt.legend()
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "01_autosomal_heterozygosity.png", dpi=300)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.scatter(valid_het["OBS_HET"], valid_het["F"],
                alpha=0.5, s=15, label="Individuals")
    if not global_outliers.empty:
        plt.scatter(global_outliers["OBS_HET"], global_outliers["F"],
                    s=25, label="Pooled heterozygosity outliers")
    plt.xlabel("Observed heterozygosity")
    plt.ylabel("PLINK inbreeding coefficient (F)")
    plt.title("Autosomal heterozygosity vs F")
    plt.legend()
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "02_heterozygosity_vs_F.png", dpi=300)
    plt.close()

    plt.figure(figsize=(9, 5))
    if not king.empty:
        plt.hist(king["KINSHIP"].dropna(), bins=100,
                 edgecolor="black", alpha=0.75)
    for threshold, label in [
        (KING_DUPLICATE, "Duplicate / MZ"),
        (KING_FIRST_DEGREE, "First degree"),
        (KING_SECOND_DEGREE, "Second degree"),
        (KING_THIRD_DEGREE, "Third degree"),
    ]:
        plt.axvline(threshold, linestyle="--",
                    label=f"{label}: {threshold}")
    plt.xlabel("KING kinship coefficient")
    plt.ylabel("Number of sample pairs")
    plt.title("KING kinship distribution")
    plt.legend(fontsize=8)
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "03_king_kinship_distribution.png", dpi=300)
    plt.close()

    plt.figure(figsize=(9, 5))
    pedsex = pd.to_numeric(sex["PEDSEX"], errors="coerce")
    for code, label in ((1, "Reported male"), (2, "Reported female")):
        vals = sex.loc[pedsex == code, "F"].dropna()
        if not vals.empty:
            plt.hist(vals, bins=40, alpha=0.55, label=label)
    plt.axvline(MAX_FEMALE_XF, linestyle="--",
                label=f"Female threshold: {MAX_FEMALE_XF}")
    plt.axvline(MIN_MALE_XF, linestyle=":",
                label=f"Male threshold: {MIN_MALE_XF}")
    plt.xlabel("Chromosome-X inbreeding coefficient (F)")
    plt.ylabel("Number of individuals")
    plt.title("Sex check: chromosome-X F-statistic")
    plt.legend()
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "04_sex_check_chrX.png", dpi=300)
    plt.close()

    plt.figure(figsize=(9, 5))
    if not y_smiss.empty:
        plt.hist(y_smiss["F_MISS"].dropna(), bins=50,
                 edgecolor="black", alpha=0.75)
    plt.axvline(Y_MISSINGNESS_THRESHOLD, linestyle="--",
                label=f"{Y_MISSINGNESS_THRESHOLD:.0%} missingness")
    plt.xlabel("Chromosome-Y missingness rate")
    plt.ylabel("Number of male individuals")
    plt.title("Chromosome-Y sample missingness")
    plt.legend()
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "05_chrY_missingness.png", dpi=300)
    plt.close()

    if population_available and "POPULATION" in valid_het.columns:
        pop_data = valid_het.dropna(subset=["POPULATION", "OBS_HET"])
        if not pop_data.empty:
            populations = sorted(pop_data["POPULATION"].unique())
            plt.figure(figsize=(max(10, len(populations) * 0.45), 6))
            values = [
                pop_data.loc[pop_data["POPULATION"] == pop, "OBS_HET"].values
                for pop in populations
            ]
            plt.boxplot(values, tick_labels=populations, showfliers=True)
            plt.xticks(rotation=90)
            plt.xlabel("1000 Genomes population")
            plt.ylabel("Observed autosomal heterozygosity")
            plt.title("Observed heterozygosity by population")
            plt.tight_layout()
            plt.savefig(PLOT_DIR / "06_heterozygosity_by_population.png",
                        dpi=300)
            plt.close()

    progress.update(1)

    add("\nPipeline completed successfully.")
    save_report()
    print("\n" + "\n".join(report))
    print(f"\nReport: {REPORT_FILE}")
    print(f"Plots:  {PLOT_DIR}")
    print(f"Logs:   {LOG_DIR}")

except Exception as exc:
    add("\nPIPELINE FAILED")
    add(repr(exc))
    save_report()
    tqdm.write(f"\nERROR: {exc}")
    tqdm.write(f"Partial report saved to: {REPORT_FILE}")
    raise

finally:
    progress.close()
