import pandas as pd


# ============================================================
# 1. LOAD DATA
# ============================================================

screening = pd.read_csv("processed/screening_long.csv")
targets = pd.read_csv("processed/targets.csv")

print("=== INPUT DATA ===")

print("\nScreening long shape:")
print(screening.shape)

print("\nScreening columns:")
print(screening.columns.tolist())

print("\nTargets shape:")
print(targets.shape)

print("\nTarget columns:")
print(targets.columns.tolist())


# ============================================================
# 2. CHECK FOR DUPLICATE COMPOUND-TARGET PAIRS
# ============================================================

duplicate_pairs = screening.duplicated(
    subset=["compound_name", "target_id"],
    keep=False
)

print("\n=== DUPLICATE CHECK ===")
print("Duplicate compound-target rows:")
print(duplicate_pairs.sum())

if duplicate_pairs.any():
    print("\nWARNING: Duplicate compound-target pairs found:")
    print(
        screening.loc[
            duplicate_pairs,
            ["compound_name", "target_id", "standardized_value"]
        ]
        .sort_values(["compound_name", "target_id"])
        .to_string(index=False)
    )


# ============================================================
# 3. MATRIX 1: FULL SCREENING MATRIX
# ============================================================

print("\n=== MATRIX 1: FULL SCREENING MATRIX ===")

screening_matrix_full = screening.pivot(
    index="compound_name",
    columns="target_id",
    values="standardized_value"
)

print("\nFull screening matrix shape:")
print(screening_matrix_full.shape)

print("\nNumber of compounds:")
print(screening_matrix_full.shape[0])

print("\nNumber of assay targets:")
print(screening_matrix_full.shape[1])

print("\nMissing cells:")
print(screening_matrix_full.isna().sum().sum())

screening_matrix_full.to_csv(
    "processed/screening_matrix_full.csv"
)

print("\nSaved:")
print("processed/screening_matrix_full.csv")


# ============================================================
# 4. INSPECT TARGET METADATA
# ============================================================

print("\n=== TARGET METADATA ===")

print("\nSpecies:")
print(targets["species"].value_counts(dropna=False))

print("\nWT / mutant classification:")
print(
    targets["wild_type_or_mutant"]
    .value_counts(dropna=False)
)

print("\nUnique biological kinase identities:")
print(targets["target_name_clean"].nunique())


# ============================================================
# 5. SELECT WT-HUMAN ASSAYS
# ============================================================

wt_human_targets = targets.loc[
    (targets["species"] == "human")
    & (targets["wild_type_or_mutant"] == "wild_type"),
    "target_id"
].tolist()

print("\n=== WT-HUMAN ASSAY SELECTION ===")

print("\nNumber of WT-human assay columns:")
print(len(wt_human_targets))


# ============================================================
# 6. FILTER SCREENING LONG DATA
# ============================================================

screening_wt_human = screening.loc[
    screening["target_id"].isin(wt_human_targets)
].copy()

print("\nWT-human long-data shape:")
print(screening_wt_human.shape)

print("\nWT-human unique compounds:")
print(screening_wt_human["compound_name"].nunique())

print("\nWT-human unique assay targets:")
print(screening_wt_human["target_id"].nunique())


# ============================================================
# 7. MATRIX 2: WT-HUMAN SCREENING MATRIX
# ============================================================

print("\n=== MATRIX 2: WT-HUMAN SCREENING MATRIX ===")

screening_matrix_wt_human = screening_wt_human.pivot(
    index="compound_name",
    columns="target_id",
    values="standardized_value"
)

print("\nWT-human screening matrix shape:")
print(screening_matrix_wt_human.shape)

print("\nNumber of compounds:")
print(screening_matrix_wt_human.shape[0])

print("\nNumber of WT-human assay columns:")
print(screening_matrix_wt_human.shape[1])

print("\nMissing cells:")
print(screening_matrix_wt_human.isna().sum().sum())


# ============================================================
# 8. CONSISTENCY CHECK
# ============================================================

expected_cells = (
    screening_matrix_wt_human.shape[0]
    * screening_matrix_wt_human.shape[1]
)

actual_rows = len(screening_wt_human)

print("\n=== MATRIX CONSISTENCY CHECK ===")

print("Expected matrix cells:")
print(expected_cells)

print("\nRows in WT-human long data:")
print(actual_rows)

if expected_cells == actual_rows:
    print("\nPASS: Matrix size matches long-data row count.")
else:
    print("\nCHECK REQUIRED: Matrix size does not match long-data row count.")


# ============================================================
# 9. SAVE MATRIX 2
# ============================================================

screening_matrix_wt_human.to_csv(
    "processed/screening_matrix_wt_human.csv"
)

print("\nSaved:")
print("processed/screening_matrix_wt_human.csv")




# ============================================================
# 10. FINAL SUMMARY
# ============================================================

print("\n=== MATRIX CREATION SUMMARY ===")

print(
    "Full screening matrix:",
    screening_matrix_full.shape
)

print(
    "WT-human screening matrix:",
    screening_matrix_wt_human.shape
)

print("\nMatrix creation complete.")

# ============================================================
# 11. LOAD Kd DATA FOR MATRIX 3
# ============================================================

print("\n=== MATRIX 3: Kd / pKd MATRIX ===")

kd = pd.read_csv("processed/kd_long.csv")

# Standardize column names for Matrix 3
kd = kd.rename(columns={
    "Entrez Gene\nSymbol": "target_id"
})

print("\nKd columns after standardization:")
print(kd.columns.tolist()) 

# Convert Kd values to numeric
kd["Kd (nM)"] = pd.to_numeric(
    kd["Kd (nM)"],
    errors="coerce"
)
print("\nKd data type after conversion:")
print(kd["Kd (nM)"].dtype)

print("\nMissing Kd after numeric conversion:")
print(kd["Kd (nM)"].isna().sum())

print("\nModifier values:")
print(kd["Modifier"].value_counts(dropna=False))


print("\nKd long shape:")
print(kd.shape)

print("\nKd columns:")
print(kd.columns.tolist())

print("\nFirst 5 Kd rows:")
print(kd.head().to_string(index=False))


# ============================================================
# 12. BASIC Kd VALIDATION
# ============================================================

print("\n=== Kd VALIDATION ===")

print("\nMissing Kd values:")
print(kd["Kd (nM)"].isna().sum())

print("\nKd values <= 0:")
print((kd["Kd (nM)"] <= 0).sum())

print("\nMinimum Kd:")
print(kd["Kd (nM)"].min())

print("\nMaximum Kd:")
print(kd["Kd (nM)"].max())

print("\nUnique compounds:")
print(kd["Compound Name"].nunique())

print("\nUnique targets:")
print(kd["target_id"].nunique())


# ============================================================
# 13. CHECK REPEATED COMPOUND-TARGET MEASUREMENTS
# ============================================================

kd_duplicates = kd.duplicated(
    subset=["Compound Name", "target_id"],
    keep=False
)

print("\n=== REPEATED Kd MEASUREMENTS ===")

print("\nNumber of rows belonging to repeated pairs:")
print(kd_duplicates.sum())

print("\nNumber of unique repeated compound-target pairs:")
print(
    kd.loc[
        kd_duplicates,
        ["Compound Name", "target_id"]
    ]
    .drop_duplicates()
    .shape[0]
)

if kd_duplicates.any():
    print("\nRepeated measurements:")
    print(
        kd.loc[
            kd_duplicates
        ]
        .sort_values(["Compound Name", "target_id"])
        .to_string(index=False)
    )

# ============================================================
# 14. INVESTIGATE REPEATED Kd PAIRS
# ============================================================

print("\n=== INVESTIGATE REPEATED Kd PAIRS ===")

repeated_detail = kd.loc[
    kd_duplicates,
    [
        "Compound Name",
        "target_id",
        "Modifier",
        "Kd (nM)"
    ]
].copy()

print("\nRepeated Kd rows:")
print(repeated_detail.to_string(index=False))


# Summarize each repeated compound-target pair
repeat_summary = (
    repeated_detail
    .groupby(["Compound Name", "target_id"])
    .agg(
        number_of_rows=("Kd (nM)", "size"),
        number_of_kd_values=("Kd (nM)", "nunique"),
        number_of_modifiers=("Modifier", "nunique")
    )
    .reset_index()
)

print("\n=== REPEATED PAIR SUMMARY ===")
print(repeat_summary.to_string(index=False))

print("\nNumber of repeated compound-target pairs:")
print(len(repeat_summary))

print("\nPairs having more than one Kd value:")
print(
    (repeat_summary["number_of_kd_values"] > 1).sum()
)

print("\nPairs having more than one modifier:")
print(
    (repeat_summary["number_of_modifiers"] > 1).sum()
)
import pandas as pd
import numpy as np
# ============================================================
# 15. CALCULATE pKd
# ============================================================

print("\n=== CALCULATE pKd ===")

# pKd = 9 - log10(Kd in nM)
kd["pkd"] = 9 - np.log10(kd["Kd (nM)"])

print("\npKd summary:")
print(kd["pkd"].describe())

print("\nMissing pKd values:")
print(kd["pkd"].isna().sum())


# ============================================================
# 16. CHECK CENSORED Kd VALUES
# ============================================================

print("\n=== Kd MODIFIERS ===")

print(kd["Modifier"].value_counts(dropna=False))
# Exact measurements have "=" or no modifier.
# Values with ">" or "<" are censored measurements.

kd["is_censored"] = kd["Modifier"].isin([">", "<"])

print("\nNumber of censored Kd measurements:")
print(kd["is_censored"].sum())

print("\nNumber of exact Kd measurements:")
print((~kd["is_censored"]).sum())

# ============================================================
# 17. AGGREGATE REPEATED EXACT MEASUREMENTS
# ============================================================

kd_for_matrix = (
    kd.groupby(
        ["Compound Name", "target_id"],
        as_index=False
    )
    .agg(
        pkd=("pkd", "mean"),
        number_of_measurements=("pkd", "size"),
        any_censored=("is_censored", "any")
    )
)

print("\nRows after resolving repeated pairs:")
print(len(kd_for_matrix))

print("\nPairs based on multiple measurements:")
print(
    (kd_for_matrix["number_of_measurements"] > 1).sum()
)
# Censored measurements should not be represented
# as exact pKd values in the matrix.
kd_for_matrix.loc[
    kd_for_matrix["any_censored"],
    "pkd"
] = np.nan
# ============================================================
# 18. MATRIX 3: pKd MATRIX
# ============================================================

print("\n=== MATRIX 3: pKd MATRIX ===")

kd_matrix_pkd = kd_for_matrix.pivot(
    index="Compound Name",
    columns="target_id",
    values="pkd"
)

print("\npKd matrix shape:")
print(kd_matrix_pkd.shape)

print("\nNumber of compounds:")
print(kd_matrix_pkd.shape[0])

print("\nNumber of targets:")
print(kd_matrix_pkd.shape[1])

print("\nNumber of measured pKd cells:")
print(kd_matrix_pkd.notna().sum().sum())

print("\nNumber of missing cells:")
print(kd_matrix_pkd.isna().sum().sum())


# ============================================================
# 19. SAVE MATRIX 3
# ============================================================

kd_matrix_pkd.to_csv(
    "processed/kd_matrix_pkd.csv"
)

print("\nSaved:")
print("processed/kd_matrix_pkd.csv")

# ============================================================
# 20. VERIFY WHY SOME COMPOUNDS ARE ABSENT FROM pKd MATRIX
# ============================================================

print("\n=== pKd MATRIX COVERAGE CHECK ===")

all_kd_compounds = set(kd["Compound Name"].dropna().unique())

matrix_compounds = set(kd_matrix_pkd.index)

missing_from_matrix = all_kd_compounds - matrix_compounds

print("\nAll compounds in Kd long data:")
print(len(all_kd_compounds))

print("\nCompounds represented in pKd matrix:")
print(len(matrix_compounds))

print("\nCompounds absent from pKd matrix:")
print(len(missing_from_matrix))

# Check whether absent compounds contain only censored measurements
missing_rows = kd[
    kd["Compound Name"].isin(missing_from_matrix)
].copy()

print("\nMeasurements belonging to absent compounds:")
print(len(missing_rows))

print("\nModifier distribution for absent compounds:")
print(missing_rows["Modifier"].value_counts(dropna=False))

print("\nAre all measurements for absent compounds censored?")
print(missing_rows["is_censored"].all())

# ============================================================
# CLEAN Kd MODIFIERS
# ============================================================

kd["modifier_clean"] = (
    kd["Modifier"]
    .astype(str)
    .str.replace("\n", "", regex=False)
    .str.replace("\r", "", regex=False)
    .str.strip()
)

print("\n=== CLEANED MODIFIER VALUES ===")
print(kd["modifier_clean"].value_counts(dropna=False))

suspicious_modifiers = ~kd["modifier_clean"].isin(
    ["=", ">", "<"]
)

print("\nNumber of suspicious modifier rows:")
print(suspicious_modifiers.sum())

print("\nSuspicious modifier rows:")
print(
    kd.loc[
        suspicious_modifiers,
        [
            "Compound Name",
            "target_id",
            "Modifier",
            "modifier_clean",
            "Kd (nM)"
        ]
    ].to_string(index=False)
)

print("\n=== ALL 18 SUSPICIOUS Kd ROWS ===")

print(
    kd.loc[
        suspicious_modifiers
    ].to_string(index=False)
)

print("\nNumber of suspicious rows:")
print(suspicious_modifiers.sum())
print("\n=== CHECK COMPONENT COMPOUND NAMES ===")

for value in kd.loc[suspicious_modifiers, "Compound Name"]:
    print("\nProblem row:", repr(value))

    parts = str(value).split("\n")

    for part in parts:
        part = part.strip()

        count = (
            kd["Compound Name"]
            .astype(str)
            .str.strip()
            .eq(part)
            .sum()
        )

        print("   ", part, "-> separate rows elsewhere:", count)
print("\n=== STRUCTURE OF 18 SUSPICIOUS ROWS ===")

suspicious_rows = kd.loc[suspicious_modifiers].copy()

for idx, row in suspicious_rows.iterrows():

    print("\n----------------------------------")
    print("ROW:", idx)

    for column in [
        "Compound Name",
        "Smiles",
        "DiscoveRx\n%Inh >90  Gene\nSymbol",
        "target_id",
        "Modifier",
        "Kd (nM)"
    ]:
        if column in kd.columns:
            value = row[column]

            print(
                column,
                "=>",
                repr(value)
            )

raise SystemExit
