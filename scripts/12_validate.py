
import pandas as pd
import numpy as np
from pathlib import Path

folder = Path("processed")

# --------------------------------------------------
# 1. LOAD DATA
# --------------------------------------------------

screening = pd.read_csv(folder / "screening_long.csv")
compounds = pd.read_csv(folder / "compounds.csv")
targets = pd.read_csv(folder / "targets.csv")

print("\n=== S4 DATASET SUMMARY ===")
print("Screening measurements:", len(screening))
print("Unique compounds:", screening["compound_name"].nunique())
print("Unique assay IDs:", screening["target_id"].nunique())

# --------------------------------------------------
# 2. MISSING VALUES AND DUPLICATES
# --------------------------------------------------

print("\n=== MISSING VALUES ===")

for column in [
    "compound_name",
    "target_id",
    "raw_value"
]:
    if column in screening.columns:
        print(column, ":", screening[column].isna().sum())

print("\n=== DUPLICATE MEASUREMENTS ===")

duplicates = screening.duplicated(
    subset=["compound_name", "target_id"],
    keep=False
)

print("Duplicate rows:", duplicates.sum())

screening[duplicates].to_csv(
    folder / "s4_duplicate_review.csv",
    index=False
)

# --------------------------------------------------
# 3. SCREENING VALUE CHECKS
# --------------------------------------------------

print("\n=== SCREENING VALUE CHECKS ===")

values = pd.to_numeric(
    screening["raw_value"],
    errors="coerce"
)

print("Non-numeric values:", values.isna().sum())
print("Minimum value:", values.min())
print("Maximum value:", values.max())

# Report values outside 0-100.
# Do not delete or clip them automatically.
outside = screening[
    values.notna() & ((values < 0) | (values > 100))
]

print("Values outside 0-100:", len(outside))

outside.to_csv(
    folder / "s4_value_review.csv",
    index=False
)

# --------------------------------------------------
# 4. TARGET CLASSIFICATION
# --------------------------------------------------

print("\n=== TARGET CLASSIFICATION ===")

print("Species:")
if "species" in targets.columns:
    print(targets["species"].value_counts(dropna=False))

print("\nWild-type / mutant:")
if "wild_type_or_mutant" in targets.columns:
    print(
        targets["wild_type_or_mutant"].value_counts(
            dropna=False
        )
    )

print("\nTarget class:")
if "target_class" in targets.columns:
    print(
        targets["target_class"].value_counts(
            dropna=False
        )
    )

# --------------------------------------------------
# 5. COMPOUND STRUCTURES
# --------------------------------------------------

print("\n=== COMPOUND STRUCTURES ===")

print("Total compound records:", len(compounds))

print("Available columns:", compounds.columns.tolist())

if "structure_status" in compounds.columns:
    print(
        compounds["structure_status"].value_counts(
            dropna=False
        )
    )

if "smiles_canonical_rdkit" in compounds.columns:
    print(
        "Missing canonical SMILES:",
        compounds["smiles_canonical_rdkit"].isna().sum()
    )

# Independently test SMILES if RDKit is installed.
try:
    from rdkit import Chem

    smiles_column = next(
        (
            col for col in [
                "smiles_original",
                "smiles_canonical_rdkit"
            ]
            if col in compounds.columns
        ),
        None
    )

    if smiles_column:
        invalid = []

        for index, smiles in compounds[
            smiles_column
        ].items():

            if pd.isna(smiles) or not str(smiles).strip():
                invalid.append(index)
            elif Chem.MolFromSmiles(str(smiles)) is None:
                invalid.append(index)

        print("Invalid or missing SMILES:", len(invalid))

        compounds.loc[invalid].to_csv(
            folder / "s4_structure_review.csv",
            index=False
        )

except ImportError:
    print("RDKit not installed: structure parsing not checked")

# --------------------------------------------------
# 6. MATRIX CONSISTENCY
# --------------------------------------------------

print("\n=== MATRIX CONSISTENCY ===")

matrix_files = [
    "screening_matrix_full.csv",
    "screening_matrix_wt_human.csv"
]

for filename in matrix_files:

    path = folder / filename

    if not path.exists():
        print(filename, ": NOT FOUND")
        continue

    matrix = pd.read_csv(path, index_col=0)

    print("\n", filename)
    print("Dimensions:", matrix.shape)
    print(
        "Non-missing cells:",
        int(matrix.notna().sum().sum())
    )

# --------------------------------------------------
# 7. SAVE SUMMARY
# --------------------------------------------------

print("\n=== S4 AUTOMATED VALIDATION FINISHED ===")
print("Review any reported problems before final approval.")
