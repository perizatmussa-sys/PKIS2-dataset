import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
# ============================================================
# REMOVE RESIDUAL MALFORMED EXTRACTION ROWS
# ============================================================

print("\n=== RESIDUAL MALFORMED ROW CHECK ===")

malformed_mask = (
    df["Kd (nM)"].isna()
    & df["Entrez Gene\nSymbol"].isna()
)

malformed_rows = df.loc[malformed_mask].copy()

print("Malformed rows found:")
print(len(malformed_rows))

print("\nMalformed compound names:")
print(
    malformed_rows["Compound Name"]
    .to_string(index=False)
)

# Remove rows that do not contain a usable target or Kd value
df = df.loc[~malformed_mask].copy()

print("\nRows after removing malformed extraction rows:")
print(len(df))

compounds = pd.read_csv("processed/compounds.csv")

print("Compound shape:")
print(compounds.columns.tolist())

print("\nFirst compounds:")
print(compounds.head())

compounds["mol"] = compounds["smiles_canonical"].apply(
    Chem.MolFromSmiles
)
print("Valid reconstruected molecules: ")
print(compounds["mol"].notna().sum())

morgan_generator = rdFingerprintGenerator.GetMorganGenerator(
    radius = 2, 
    fpSize = 2048
)
first_mol = compounds["mol"].iloc[0]

first_fp = morgan_generator.GetFingerprint(first_mol)

print("\nFirst fingerprint:")
print(first_fp)
# ============================================================
# FINAL CLEANUP OF Kd DATA
# ============================================================

print("\n=== FINAL Kd CLEANUP ===")

print("\nRows before final cleanup:")
print(len(df))


# ------------------------------------------------------------
# 1. Find residual malformed extraction rows
# ------------------------------------------------------------

malformed_mask = (
    df["Kd (nM)"].isna()
    & df["Entrez Gene\nSymbol"].isna()
)

malformed_rows = df.loc[malformed_mask].copy()

print("\nMalformed extraction rows found:")
print(len(malformed_rows))

if len(malformed_rows) > 0:
    print("\nMalformed rows:")
    print(
        malformed_rows.to_string(index=False)
    )


# ------------------------------------------------------------
# 2. SAFETY CHECK
# ------------------------------------------------------------

# We currently expect exactly 18 malformed extraction rows.
# Stop the script if a different number is found.
expected_malformed_rows = 18

if len(malformed_rows) != expected_malformed_rows:
    raise ValueError(
        f"Expected {expected_malformed_rows} malformed rows, "
        f"but found {len(malformed_rows)}. "
        "STOPPING so that valid Kd observations are not removed."
    )

print(
    "\nPASS: Expected number of malformed extraction rows found."
)


# ------------------------------------------------------------
# 3. Remove ONLY those malformed rows
# ------------------------------------------------------------

df = df.loc[~malformed_mask].copy()

df = df.reset_index(drop=True)

print("\nRows after removing malformed extraction rows:")
print(len(df))


# ------------------------------------------------------------
# 4. Convert Kd to numeric
# ------------------------------------------------------------

df["Kd (nM)"] = pd.to_numeric(
    df["Kd (nM)"],
    errors="coerce"
)


# ------------------------------------------------------------
# 5. Final missing-value checks
# ------------------------------------------------------------

print("\n=== FINAL MISSING-VALUE CHECK ===")

print("\nMissing compound names:")
print(df["Compound Name"].isna().sum())

print("\nMissing target symbols:")
print(df["Entrez Gene\nSymbol"].isna().sum())

print("\nMissing Kd values:")
print(df["Kd (nM)"].isna().sum())


# ------------------------------------------------------------
# 6. Kd validity checks
# ------------------------------------------------------------

print("\n=== Kd VALUE CHECK ===")

print("\nKd values <= 0:")
print((df["Kd (nM)"] <= 0).sum())

print("\nMinimum Kd:")
print(df["Kd (nM)"].min())

print("\nMaximum Kd:")
print(df["Kd (nM)"].max())


# ------------------------------------------------------------
# 7. Modifier inspection
# ------------------------------------------------------------

print("\n=== MODIFIER CHECK ===")

print(df["Modifier"].value_counts(dropna=False))


# ------------------------------------------------------------
# 8. Dataset summary
# ------------------------------------------------------------

print("\n=== FINAL Kd DATASET SUMMARY ===")

print("Rows:")
print(len(df))

print("\nUnique compounds:")
print(df["Compound Name"].nunique())

print("\nUnique targets:")
print(df["Entrez Gene\nSymbol"].nunique())


# ------------------------------------------------------------
# 9. Save cleaned Kd long dataset
# ------------------------------------------------------------

df.to_csv(
    "processed/kd_long.csv",
    index=False
)

print("\nSaved cleaned Kd dataset:")
print("processed/kd_long.csv")