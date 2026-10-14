import pandas as pd

file_path = "raw/pone.0181585.s003.xlsx"

excel_file = pd.ExcelFile(file_path)

print("\nSheet names:")
print(excel_file.sheet_names)

all_tables = []

df = pd.read_excel(
    file_path,
    sheet_name="Table 1"
)

column_names = df.columns

all_tables.append(df)

for sheet in excel_file.sheet_names[1:]:

    temp_df = pd.read_excel(
        file_path,
        sheet_name=sheet,
        header=None
    )

    temp_df.columns = column_names

    all_tables.append(temp_df)

combined_df = pd.concat(
    all_tables,
    ignore_index=True
)

print("\nALL TABLES COMBINED:")
print(combined_df.shape)

print("\nFirst 10 rows of combined data:")
print(combined_df.head(10))

print("\nLast 10 rows of combined data:")
print(combined_df.tail(10))

print("\nMissing values in combined data:")
print(combined_df.isna().sum())

print("\nCombined columns:")
print(combined_df.columns.tolist())

missing_compound = combined_df[
    combined_df["Compound Name"].isna()
]
print("\nRows with missing Compound Name:")
print(missing_compound.head(20))

print("\nRows 0 to 25:")
print(
    combined_df.loc[
        0:25,
        ["Compound Name", "Smiles", "Entrez Gene\nSymbol", "Modifier", "Kd (nM)"]
    ].to_string()
)
clean_df = combined_df.copy()

continuation_rows = (
    clean_df["Compound Name"].isna()
    & clean_df["Kd (nM)"].isna()
    & clean_df["Smiles"].notna()
)
print("\nPossible continuation rows:")
print(continuation_rows.sum())

continuation_indices = clean_df.index[continuation_rows]

print("\nFirst 20 continuation indices:")
print(continuation_indices[:20])

for i in continuation_indices[:20]:
    print(
        i,
        "-> next row:",
        i+1, 
        "| Compound:",
        clean_df.loc[i+1, "Compound Name"]
    )

    bad_continuations = []

for i in continuation_indices:

    if i + 1 >= len(clean_df):
        bad_continuations.append(i)

    elif pd.isna(clean_df.loc[i + 1, "Compound Name"]):
        bad_continuations.append(i)

print("\nContinuation rows without a valid next compound:")
print(len(bad_continuations))

print("\nTheir row numbers:")
print(bad_continuations)

for i in bad_continuations:
    print("\n--------------------")
    print("Problem around row:", i)

    start = max(0, i - 2)
    end = min(len(clean_df) - 1, i + 3)

    print(
        clean_df.loc[
            start:end,
            ["Compound Name", "Smiles", "Entrez Gene\nSymbol", "Modifier", "Kd (nM)"]
        ].to_string()
    )
    max_run = 0
current_run = 0

for is_continuation in continuation_rows:

    if is_continuation:
        current_run += 1
        max_run = max(max_run, current_run)
    else:
        current_run = 0

print("\nMaximum consecutive continuation rows:")
print(max_run)

smiles_parts = []
rows_to_remove = []

for i, row in clean_df.iterrows():
    if pd.isna(row["Compound Name"]):
                smiles_parts.append(row["Smiles"])
                rows_to_remove.append(i)
    else:
       if smiles_parts:
            current_smiles = row["Smiles"]

            if pd.notna(current_smiles):
                 full_smiles = "".join(str(part) for part in smiles_parts) + str(current_smiles)
                 clean_df.at[i, "Smiles"] = full_smiles
       smiles_parts = []
print("\nNumber of continuation rows marked for removal:")
print(len(rows_to_remove))

clean_df = clean_df.drop(index=rows_to_remove)
clean_df = clean_df.reset_index(drop=True)
print("\nShape after removing continuation rows:")
print(clean_df.shape)
print("\nMissing values after cleaning:")
print(clean_df.isna().sum())
print("\nFirst 10 cleaned rows:")
print(clean_df.head(10))
print("\nMissing Compound Names after cleaning:")
print(clean_df["Compound Name"].isna().sum())
print("\nMissing Kd values after cleaning:")
print(clean_df["Kd (nM)"].isna().sum())
print("\nFinal shape:")
print(clean_df.shape)

missing_kd = clean_df[clean_df["Kd (nM)"].isna()]

print("\nRows with missing Kd")
print(missing_kd.to_string())

missing_kd_indices = clean_df.index[clean_df["Kd (nM)"].isna()]

for i in missing_kd_indices:
     print("\n-----------------")
     print("Missing Kd around row:", i)

     start = max(0, i -1 )
     end = min(len(clean_df) - 1, i + 1)

     print(
          clean_df.loc[
               start:end,
               ["Compound Name", "Smiles", "Entrez Gene\nSymbol", "Modifier", "Kd (nM)" ]
          ].to_string()
     )
print("\nSpecial missing-Kd rows:")
print(
    missing_kd[
        ["Compound Name", "Smiles",
         "Entrez Gene\nSymbol", "Modifier", "Kd (nM)"]
    ].to_string()
)
special_rows = (
    clean_df["Kd (nM)"].isna()
    & clean_df["Smiles"].isna()
    & clean_df["Entrez Gene\nSymbol"].isna()
    & clean_df["Modifier"].isna()
)

print("\nNumber of special artifact rows:")
print(special_rows.sum())
clean_df = clean_df[~special_rows].copy()
clean_df = clean_df.reset_index(drop=True)

print("\nFinal shape:")
print(clean_df.shape)

print("\nMissing Compound Names:")
print(clean_df["Compound Name"].isna().sum())

print("\nMissing Kd values:")
print(clean_df["Kd (nM)"].isna().sum())

output_path = "processed/kd_long.csv"

clean_df.to_csv(
     output_path,
     index = False
)
print("\nSaved cleaned Kd data to:")
print(output_path)