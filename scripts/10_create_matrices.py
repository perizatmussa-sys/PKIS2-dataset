
import pandas as pd
import numpy as np
from pathlib import Path

folder = Path("processed")

# Load the latest candidate Kd measurements.
kd = pd.read_csv(folder / "kd_long_candidate.csv")

print("=== INPUT DATA ===")
print("Candidate Kd rows:", len(kd))

# Only exact Kd measurements can be converted
# into exact pKd values.
exact = kd[kd["kd_operator"] == "="].copy()

exact["kd_nm"] = pd.to_numeric(
    exact["kd_nm"], errors="coerce"
)

exact = exact[
    exact["kd_nm"].notna()
    & (exact["kd_nm"] > 0)
].copy()

exact["pkd"] = 9 - np.log10(exact["kd_nm"])

# Check duplicate compound-target pairs.
duplicates = exact[
    exact.duplicated(
        subset=["compound_name", "target_id"],
        keep=False
    )
].copy()

duplicates.to_csv(
    folder / "kd_matrix_duplicate_review.csv",
    index=False
)

# Exclude ambiguous duplicate pairs for now.
# Do not average measurements without verification.
ambiguous_pairs = duplicates[
    ["compound_name", "target_id"]
].drop_duplicates()

if not ambiguous_pairs.empty:
    exact = exact.merge(
        ambiguous_pairs.assign(ambiguous=True),
        on=["compound_name", "target_id"],
        how="left"
    )

    exact = exact[
        exact["ambiguous"].isna()
    ].drop(columns="ambiguous")

# Create compound x kinase matrix.
matrix = exact.pivot(
    index="compound_name",
    columns="target_id",
    values="pkd"
)

matrix.to_csv(
    folder / "kd_matrix_pkd_candidate.csv"
)

print("\n=== CANDIDATE KD MATRIX ===")
print("Exact Kd measurements used:", len(exact))
print("Compounds:", matrix.shape[0])
print("Targets:", matrix.shape[1])
print("Duplicate rows excluded:", len(duplicates))

print("\nFiles saved:")
print("processed/kd_matrix_pkd_candidate.csv")
print("processed/kd_matrix_duplicate_review.csv")
