
import pandas as pd
from pathlib import Path

folder = Path("processed")

s3 = pd.read_csv(folder / "kd_long_candidate.csv", dtype=str)
s4 = pd.read_csv(folder / "screening_long.csv", dtype=str)

# Standardize identifiers without modifying original files.
for df in [s3, s4]:
    df["compound_name"] = df["compound_name"].str.strip().str.upper()
    df["target_id"] = df["target_id"].str.strip().str.upper()

# Unique compound and target identifiers in each dataset.
s3_compounds = set(s3["compound_name"].dropna())
s4_compounds = set(s4["compound_name"].dropna())

s3_targets = set(s3["target_id"].dropna())
s4_targets = set(s4["target_id"].dropna())

matched_compounds = s3_compounds & s4_compounds
matched_targets = s3_targets & s4_targets

unresolved_compounds = sorted(s3_compounds - s4_compounds)
unresolved_targets = sorted(s3_targets - s4_targets)

# Check repeated compound-target combinations.
duplicates = s3[
    s3.duplicated(
        subset=["compound_name", "target_id"],
        keep=False
    )
].copy()

# Save reports.
pd.DataFrame({
    "compound_name": unresolved_compounds
}).to_csv(folder / "s3_unresolved_compounds.csv", index=False)

pd.DataFrame({
    "target_id": unresolved_targets
}).to_csv(folder / "s3_unresolved_targets.csv", index=False)

duplicates.to_csv(folder / "s3_duplicate_pairs.csv", index=False)

# Check compound-target pair overlap without combining measurements.
s3_pairs = set(zip(s3["compound_name"], s3["target_id"]))
s4_pairs = set(zip(s4["compound_name"], s4["target_id"]))

print("\n=== S3-S4 MAPPING REPORT ===")

print("S3 unique compounds:", len(s3_compounds))
print("S4 unique compounds:", len(s4_compounds))
print("Matched compound identifiers:", len(matched_compounds))
print("Unresolved compound identifiers:", len(unresolved_compounds))

print("\nS3 unique targets:", len(s3_targets))
print("S4 unique targets:", len(s4_targets))
print("Matched target identifiers:", len(matched_targets))
print("Unresolved target identifiers:", len(unresolved_targets))

print("\nMatched compound-target pairs:", len(s3_pairs & s4_pairs))
print("Rows with duplicate S3 pairs:", len(duplicates))

print("\nReports saved in processed/")
