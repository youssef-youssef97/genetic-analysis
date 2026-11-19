import pandas as pd

path = "snp_array/results/QC/SexCheck/sex_check_thresholds.sexcheck"
df = pd.read_csv(path, sep=r"\s+")

print("=== STATUS SUMMARY ===")
print(df["STATUS"].value_counts())

problems = df[df["STATUS"] != "OK"].copy()

print("\n=== ALL PROBLEM SAMPLES ===")
print(problems.to_string(index=False))

problems.to_csv(
    "snp_array/results/QC/SexCheck/problem_samples.csv",
    index=False
)

print("\nSaved problem_samples.csv")