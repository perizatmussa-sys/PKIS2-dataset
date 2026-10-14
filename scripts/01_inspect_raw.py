import pandas as pd

file_path = "raw/pone.0181585.s004.xlsx"

excel_file = pd.ExcelFile(file_path)

print(excel_file.sheet_names)
df = pd.read_excel(file_path, sheet_name=excel_file.sheet_names[0])
print(df.shape)
print("\Column names:")
print(df.columns.tolist())
print("\nFirst 5 rows:")
print(df.head())
print(df.iloc[:5, :8])
print("\nLast 5 rows: ")
print(df.iloc[-5:, :8])
print("\nCompletely empty rows:")
print(df[df.isna().all(axis=1)])
print("\nNon-empty values in last row: ")
print(df.iloc[-1].dropna())
print("\n\n===== INSPECTING ALL S4 SHEETS =====")

for sheet in excel_file.sheet_names:
    print("\n==========================")
    print("SHEET:", sheet)
    print("==========================")

    temp_df = pd.read_excel(
        file_path,
        sheet_name=sheet,
        header=None
    )

    print("Shape:")
    print(temp_df.shape)

    print("\nFirst 15 rows:")
    print(temp_df.head(15).to_string())

