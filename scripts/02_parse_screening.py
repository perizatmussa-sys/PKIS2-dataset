import pandas as pd
from pathlib import Path


# ============================================================
# 1. FILE PATHS
# ============================================================

file_path = "raw/pone.0181585.s004.xlsx"
output_path = "processed/screening_long.csv"

Path("processed").mkdir(exist_ok=True)


# ============================================================
# 2. OPEN S4
# ============================================================

excel_file = pd.ExcelFile(file_path)

print("\n========================================")
print("S4 SCREENING DATA")
print("========================================")

print("\nSheet names:")
print(excel_file.sheet_names)


# ============================================================
# 3. PROCESS EACH SHEET
# ============================================================

all_long = []

for sheet in excel_file.sheet_names:

    print("\n----------------------------------------")
    print("Processing sheet:", sheet)
    print("----------------------------------------")

    df = pd.read_excel(
        file_path,
        sheet_name=sheet
    )

    print("\nOriginal shape:")
    print(df.shape)

    print("\nFirst 5 rows:")
    print(df.head().to_string())


    # ========================================================
    # 4. IDENTIFY METADATA COLUMNS
    # ========================================================

    metadata_candidates = [
        "Regno",
        "Compound",
        "Chemotype",
        "Smile"
    ]

    metadata_columns = [
        col for col in metadata_candidates
        if col in df.columns
    ]

    print("\nMetadata columns:")
    print(metadata_columns)


    # ========================================================
    # 5. IDENTIFY KINASE/TARGET COLUMNS
    # ========================================================

    target_columns = [
        col for col in df.columns
        if col not in metadata_columns
    ]

    print("\nNumber of target columns:")
    print(len(target_columns))

    print("\nFirst 20 target columns:")
    print(target_columns[:20])


    # ========================================================
    # 6. CONVERT WIDE DATA TO LONG DATA
    # ========================================================

    long_df = df.melt(
        id_vars=metadata_columns,
        value_vars=target_columns,
        var_name="target_original",
        value_name="raw_value"
    )

    print("\nShape after converting to long format:")
    print(long_df.shape)

    print("\nFirst 10 long-format rows:")
    print(long_df.head(10).to_string(index=False))


    # ========================================================
    # 7. REMOVE ROWS WITHOUT MEASUREMENTS
    # ========================================================

    print("\nMissing raw values before cleaning:")
    print(long_df["raw_value"].isna().sum())

    long_df = long_df[
        long_df["raw_value"].notna()
    ].copy()

    print("\nShape after removing missing measurements:")
    print(long_df.shape)


    # ========================================================
    # 8. MAKE RAW VALUE NUMERIC
    # ========================================================

    long_df["raw_value"] = pd.to_numeric(
        long_df["raw_value"],
        errors="coerce"
    )

    print("\nValues that could not be converted to numeric:")
    print(long_df["raw_value"].isna().sum())

    # Remove entries that were not actual numerical measurements
    long_df = long_df[
        long_df["raw_value"].notna()
    ].copy()


    # ========================================================
    # 9. STANDARDIZE BASIC COLUMN NAMES
    # ========================================================

    rename_map = {}

    if "Regno" in long_df.columns:
        rename_map["Regno"] = "compound_id"

    if "Compound" in long_df.columns:
        rename_map["Compound"] = "compound_name"

    if "Smile" in long_df.columns:
        rename_map["Smile"] = "smiles"

    if "Chemotype" in long_df.columns:
        rename_map["Chemotype"] = "chemotype"

    long_df = long_df.rename(columns=rename_map)


    # ========================================================
    # 10. CLEAN TARGET NAMES
    # ========================================================

    long_df["target_original"] = (
        long_df["target_original"]
        .astype(str)
        .str.strip()
    )

    # We are NOT mapping aliases yet.
    # For now target_id preserves the original target name.
    long_df["target_id"] = long_df["target_original"]


    # ========================================================
    # 11. PRESERVE MEASUREMENT SEMANTICS
    # ========================================================

    # IMPORTANT:
    # Do not automatically assume these measurements are
    # % inhibition or % control.
    #
    # Until the supplementary material is verified, preserve
    # the original measurement without mathematical conversion.

    long_df["measurement_type_raw"] = "screening_value"

    long_df["measurement_type_standard"] = pd.NA

    long_df["value_standard"] = pd.NA

    long_df["screening_concentration_nM"] = pd.NA

    long_df["target_class"] = "kinase"


    # ========================================================
    # 12. ADD PROVENANCE
    # ========================================================

    long_df["source_file"] = "pone.0181585.s004.xlsx"

    long_df["source_sheet"] = sheet


    # ========================================================
    # 13. STORE THIS SHEET
    # ========================================================

    all_long.append(long_df)


# ============================================================
# 14. COMBINE ALL SHEETS
# ============================================================

screening_long = pd.concat(
    all_long,
    ignore_index=True
)

print("\n========================================")
print("COMBINED SCREENING DATA")
print("========================================")

print("\nFinal shape:")
print(screening_long.shape)


# ============================================================
# 15. CHECK IMPORTANT COLUMNS
# ============================================================

print("\nMissing values:")

for column in [
    "compound_id",
    "compound_name",
    "target_id",
    "target_original",
    "raw_value"
]:
    if column in screening_long.columns:
        print(
            column,
            ":",
            screening_long[column].isna().sum()
        )


# ============================================================
# 16. CHECK DUPLICATES
# ============================================================

duplicate_count = screening_long.duplicated(
    subset=[
        "compound_id",
        "target_original"
    ]
).sum()

print("\nDuplicate compound-target pairs:")
print(duplicate_count)


# ============================================================
# 17. INSPECT MEASUREMENT DISTRIBUTION
# ============================================================

print("\nRaw value summary:")
print(
    screening_long["raw_value"].describe()
)

print("\nMinimum raw value:")
print(screening_long["raw_value"].min())

print("\nMaximum raw value:")
print(screening_long["raw_value"].max())


# ============================================================
# 18. SELECT FINAL COLUMN ORDER
# ============================================================

final_columns = [
    "compound_id",
    "compound_name",
    "target_id",
    "target_original",
    "raw_value",
    "measurement_type_raw",
    "measurement_type_standard",
    "value_standard",
    "screening_concentration_nM",
    "target_class",
    "source_file",
    "source_sheet"
]

# Keep only columns that exist
final_columns = [
    col for col in final_columns
    if col in screening_long.columns
]

screening_long = screening_long[
    final_columns
].copy()


# ============================================================
# 19. FINAL CHECK
# ============================================================

print("\n========================================")
print("FINAL SCREENING DATA")
print("========================================")

print("\nFinal shape:")
print(screening_long.shape)

print("\nFinal columns:")
print(screening_long.columns.tolist())

print("\nFirst 20 rows:")
print(
    screening_long
    .head(20)
    .to_string(index=False)
)


# ============================================================
# 20. SAVE
# ============================================================

screening_long.to_csv(
    output_path,
    index=False
)

print("\n========================================")
print("SAVED SUCCESSFULLY")
print("========================================")

print("\nOutput file:")
print(output_path)

print("\nNumber of rows saved:")
print(len(screening_long))