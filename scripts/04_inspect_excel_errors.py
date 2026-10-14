import pandas as pd

file_path = "raw/pone.0181585.s003.xlsx"

excel = pd.ExcelFile(file_path)
tables = []

for sheet_number, sheet_name in enumerate(excel.sheet_names):

    if sheet_number == 0:
        df = pd.read_excel(file_path, sheet_name=sheet_name)
        columns = df.columns.tolist()
    else:
        df = pd.read_excel(
            file_path,
            sheet_name=sheet_name,
            header=None
        )
        df.columns = columns

    df["source_sheet"] = sheet_name
    df["source_row"] = df.index + (2 if sheet_number == 0 else 1)

    tables.append(df)

data = pd.concat(tables, ignore_index=True)

# Detect suspicious records
modifier = data["Modifier"].fillna("").astype(str)
compound = data["Compound Name"].fillna("").astype(str)
kd_value = data["Kd (nM)"].fillna("").astype(str)

problems = data[
    data["Compound Name"].isna()
    | data["Entrez Gene\nSymbol"].isna()
    | data["Kd (nM)"].isna()
    | modifier.str.contains("\n", regex=False)
    | compound.str.contains("\n", regex=False)
    | kd_value.str.contains("\n", regex=False)
].copy()

problems.to_csv(
    "processed/kd_excel_conversion_review.csv",
    index=False
)

print("Total Excel rows:", len(data))
print("Rows needing review:", len(problems))
print("Missing compound names:", data["Compound Name"].isna().sum())
print("Missing target symbols:", data["Entrez Gene\nSymbol"].isna().sum())
print("Missing Kd values:", data["Kd (nM)"].isna().sum())
print("Merged modifier cells:", modifier.str.contains("\n", regex=False).sum())
print("Merged Kd cells:", kd_value.str.contains("\n", regex=False).sum())

print("\nSaved: processed/kd_excel_conversion_review.csv")
print("Original kd_long.csv was NOT changed.")

missing_compound = data["Compound Name"].isna()

pure_continuation = (
    missing_compound 
    & data["Smiles"].notna()
    & data["Entrez Gene\nSymbol"].isna()
    & data["Modifier"].isna()
    & data["Kd (nM)"].isna()
)
suspicious_continuation = (
    missing_compound
    & (
        data["Entrez Gene\nSymbol"].notna()
        | data["Modifier"].notna()
        | data["Kd (nM)"].notna()
    )
)

missing_target_measurements = (
    data["Compound Name"].notna()
    & data["Entrez Gene\nSymbol"].isna()
)
print("Likely SMILES continuation rows:", pure_continuation.sum())
print("Suspicious continuation rows:", suspicious_continuation.sum())
print("Compound rows missing target symbols:", missing_target_measurements.sum())

# Save suspicious rows for further inspection
priority_review = data.loc[
    suspicious_continuation | missing_target_measurements
].copy()

priority_review.to_csv(
    "processed/kd_priority_review.csv",
    index=False
)

print("\nSaved: processed/kd_priority_review.csv")
print("No original files were modified.")
print("\n=== INSPECT MISSING TARGET ROWS ===")

missing_target_rows = data.loc[missing_target_measurements]

print("Total missing target rows:", len(missing_target_rows))

print(
    "Missing target WITH Kd:",
    missing_target_rows["Kd (nM)"].notna().sum()
)

print(
    "Missing target WITHOUT Kd:",
    missing_target_rows["Kd (nM)"].isna().sum()
)

print(
    "Missing target WITH modifier:",
    missing_target_rows["Modifier"].notna().sum()
)

print(
    "Missing target WITH SMILES:",
    missing_target_rows["Smiles"].notna().sum()
)

print("\n=== FIRST 10 MISSING TARGET ROWS ===")

columns_to_show = [
    "Compound Name",
    "Entrez Gene\nSymbol",
    "Modifier",
    "Kd (nM)",
    "source_sheet",
    "source_row"
]

print(
    missing_target_rows[columns_to_show]
    .head(10)
    .to_string(index=False)
)

print("\n=== ALL COLUMN NAMES ===")
for column in data.columns:
    print(repr(column))

print("\n=== ALL COLUMN NAMES ===")

for column in data.columns:
    print(repr(column))


print("\n=== COMPARE GENE COLUMNS ===")

gene_columns = [
    "Compound Name",
    "DiscoverX\n%Inh >90  Gene\nSymbol",
    "Entrez Gene\nSymbol",
    "Kd (nM)",
    "source_sheet",
    "source_row"
]

print(
    missing_target_rows[gene_columns]
    .head(15)
    .to_string(index=False)
)