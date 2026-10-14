
import pandas as pd
import numpy as np
from pathlib import Path

folder = Path("processed")

kd = pd.read_csv(folder / "kd_long_candidate.csv")
matrix = pd.read_csv(
    folder / "kd_matrix_pkd_candidate.csv",
    index_col=0
)

print("=== KD VALIDATION ===")

# 1. Check missing identifiers.
print("Missing compound names:",
      kd["compound_name"].isna().sum())

print("Missing target IDs:",
      kd["target_id"].isna().sum())

# 2. Check Kd values.
kd["kd_nm"] = pd.to_numeric(kd["kd_nm"], errors="coerce")

print("Invalid or missing Kd values:",
      (kd["kd_nm"].isna() | (kd["kd_nm"] <= 0)).sum())

# 3. Check pKd calculations for exact measurements.
exact = kd[kd["kd_operator"] == "="].copy()

expected_pkd = 9 - np.log10(exact["kd_nm"])

differences = np.abs(
    pd.to_numeric(exact["pkd"], errors="coerce")
    - expected_pkd
)

print("Incorrect pKd calculations:",
      (~np.isclose(differences, 0, atol=1e-6)).sum())

# 4. Check duplicates.
print("Duplicate compound-target rows:",
      kd.duplicated(
          ["compound_name", "target_id"],
          keep=False
      ).sum())

# 5. Check matrix values.
print("Matrix dimensions:", matrix.shape)
print("Non-missing matrix values:",
      matrix.notna().sum().sum())

print("\n=== VALIDATION COMPLETE ===")
