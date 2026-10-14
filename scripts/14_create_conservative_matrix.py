
import pandas as pd
from pathlib import Path

folder = Path("processed")

full = pd.read_csv(
    folder / "screening_matrix_full.csv",
    index_col=0
)

special = pd.read_csv(
    folder / "s4_special_assays_review.csv"
)

# Use the exact S4 assay identifiers.
exclude = set(special["target_id"].dropna().astype(str))

# Keep all other assays.
conservative = full.loc[
    :,
    ~full.columns.isin(exclude)
]

conservative.to_csv(
    folder / "screening_matrix_conservative_candidate.csv"
)

print("\n=== CONSERVATIVE MATRIX ===")
print("Original assay columns:", full.shape[1])
print("Flagged assay identifiers:", len(exclude))
print("Actually excluded:", full.shape[1] - conservative.shape[1])
print("Remaining assay columns:", conservative.shape[1])
print("Compounds:", conservative.shape[0])

print("\nSaved:")
print("processed/screening_matrix_conservative_candidate.csv")
