import pandas as pd 
kd = pd.read_csv("processed/kd_long.csv")

print("\n=== Kd CLEANUP INSPECTION ===")

print("Original rows:", len(kd))

modifier = kd["Modifier"].fillna("").astype(str)

malformed_mask = modifier.str.contains(
    r"[\r\n]",
    regex=True
)
malformed_rows = kd.loc[malformed_mask].copy()
print("\nMalformed rows found:")
print(len(malformed_rows))

print("\nMissing targets among malformed rows:")
print(malformed_rows["Entrez Gene\nSymbol"].isna().sum())

print("\nMissing Kd values among malformed rows:")
print(malformed_rows["Kd (nM)"].isna().sum())

print("\nMalformed rows:")
print(
    malformed_rows[
        ["Compound Name", "Modifier", "Kd (nM)"]
    ].to_string(index=False)
)

# Stop here: inspection only
print("\nInspection complete. No files were changed.")
# ============================================================
# SAVE MALFORMED ROWS FOR MANUAL VERIFICATION
# ============================================================

malformed_rows.to_csv(
    "processed/kd_malformed_review.csv",
    index=False
)

print("\nSaved malformed rows for review:")
print("processed/kd_malformed_review.csv")

print("\nNumber of rows requiring review:")
print(len(malformed_rows))