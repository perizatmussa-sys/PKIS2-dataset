import pandas as pd

# Read the CSV extracted from the PDF
df = pd.read_csv(
    "processed/kd_pdf_raw_extracted.csv",
    header=None,
    dtype=str
)

print("=== DATASET SHAPE ===")
print("Rows:", len(df))
print("Columns:", len(df.columns))

print("\n=== FIRST 15 ROWS ===")
print(df.head(15).to_string(index=True))

print("\n=== NON-MISSING VALUES PER COLUMN ===")
print(df.notna().sum())

print("\n=== COLUMN VALUE EXAMPLES ===")

for column in df.columns:
    print(f"\nColumn {column}:")

    examples = df[column].dropna().head(5)

    for value in examples:
        print(repr(value[:100]))

print("\n=== ROWS WITH EXTRA COLUMNS ===")

for column in [5, 6]:
    print(f"\nColumn {column}")

    rows = df[df[column].notna()]

    print("Total rows:", len(rows))

    for index, row in rows.head(10).iterrows():
        print(f"\nRow {index}")

        for col in df.columns:
            value = row[col]

            if pd.notna(value):
                print(f"Column {col}: {str(value)[:100]}")

print("\n=== PDF LAYOUT COUNTS ===")

layout_a = df[
    df[4].notna()
    & df[5].isna()
    & df[6].isna()
]

layout_b = df[
    df[5].notna()
    & df[6].notna()
]

print("Layout A rows:", len(layout_a))
print("Layout B rows:", len(layout_b))

print("\n=== KD OPERATORS IN LAYOUT B ===")
print(layout_b[5].value_counts(dropna=False))