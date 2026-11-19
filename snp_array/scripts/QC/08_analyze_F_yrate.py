import pandas as pd

path = "snp_array/results/QC/SexCheck/sex_check.sexcheck"

df = pd.read_csv(path, sep=r"\s+")

for sex, label in [(1, "Males"), (2, "Females")]:
    subset = df[df["PEDSEX"] == sex]

    print(f"\n{'='*15} {label} {'='*15}")
    print(subset[["F", "YRATE"]].describe(
        percentiles=[0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99]
    ).round(4))

print("\nMissing values:")
print(df[["F", "YRATE"]].isna().sum())