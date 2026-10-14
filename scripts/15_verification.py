
import pandas as pd
from pathlib import Path

# ----------------------------------------
# 1. FILE PATHS
# ----------------------------------------

folder = Path("processed")

s3 = pd.read_csv(folder / "kd_long_candidate.csv")
s4 = pd.read_csv(folder / "screening_long.csv")

# Fixed seed ensures reproducible sampling.
SEED = 42
SAMPLE_SIZE = 20

# ----------------------------------------
# 2. S3 VERIFICATION SAMPLE
# ----------------------------------------

# Select 20 candidate Kd measurements.
s3_sample = s3.sample(
    n=min(SAMPLE_SIZE, len(s3)),
    random_state=SEED
).copy()

# Keep source row to help locate PDF records.
s3_columns = [
    "source_row",
    "compound_name",
    "target_original",
    "target_id",
    "kd_operator",
    "kd_nm",
    "pkd"
]

s3_sample = s3_sample[
    [col for col in s3_columns if col in s3_sample.columns]
]

s3_sample["source_value_checked"] = ""
s3_sample["verification_status"] = "not_checked"
s3_sample["verification_notes"] = ""

s3_sample.to_csv(
    folder / "s3_manual_verification_sample.csv",
    index=False
)

# ----------------------------------------
# 3. S4 VERIFICATION SAMPLE
# ----------------------------------------

# Select 20 screening measurements.
s4_sample = s4.sample(
    n=min(SAMPLE_SIZE, len(s4)),
    random_state=SEED
).copy()

s4_columns = [
    "compound_name",
    "target_id",
    "target_original",
    "raw_value",
    "raw_measurement_type",
    "standardized_value",
    "standardized_measurement_type",
    "source_file",
    "source_sheet"
]

s4_sample = s4_sample[
    [col for col in s4_columns if col in s4_sample.columns]
]

s4_sample["source_value_checked"] = ""
s4_sample["verification_status"] = "not_checked"
s4_sample["verification_notes"] = ""

s4_sample.to_csv(
    folder / "s4_manual_verification_sample.csv",
    index=False
)

# ----------------------------------------
# 4. PRINT SUMMARY
# ----------------------------------------

print("\n=== MANUAL VERIFICATION SAMPLES ===")
print("S3 measurements selected:", len(s3_sample))
print("S4 measurements selected:", len(s4_sample))

print("\n=== S3 SAMPLE PREVIEW ===")
print(
    s3_sample[
        ["compound_name", "target_id", "kd_nm"]
    ].head(5).to_string(index=False)
)

print("\n=== S4 SAMPLE PREVIEW ===")
print(
    s4_sample[
        ["compound_name", "target_id", "raw_value"]
    ].head(5).to_string(index=False)
)

print("\nFiles saved:")
print("processed/s3_manual_verification_sample.csv")
print("processed/s4_manual_verification_sample.csv")
