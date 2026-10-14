import pdfplumber

pdf_path = "raw/pone.0181585.s003.pdf"

table_settings = {
    "vertical_strategy": "text",
    "horizontal_strategy": "text",
    "intersection_tolerance": 5,
}

all_rows = []

with pdfplumber.open(pdf_path) as pdf:

    print("Total PDF pages:", len(pdf.pages))

    for page_number, page in enumerate(pdf.pages, start=1):

        tables = page.extract_tables(table_settings)

        page_rows = 0

        for table in tables:
            for row in table:
                all_rows.append(row)
                page_rows += 1

        print(
            f"Page {page_number}: "
            f"{len(tables)} tables, "
            f"{page_rows} rows"
        )

print("\n=== EXTRACTION SUMMARY ===")

print("Total extracted rows:", len(all_rows))

print("\nFirst 5 extracted rows:")

for row in all_rows[:5]:
    print(row)

import pandas as pd

# Find the maximum number of columns
max_columns = max(len(row) for row in all_rows)

# Make every row the same length
normalized_rows = []

for row in all_rows:
    new_row = list(row)
    while len(new_row) < max_columns:
        new_row.append(None)
    normalized_rows.append(new_row)

# Create a DataFrame
pdf_df = pd.DataFrame(normalized_rows)

# Save extracted data
pdf_df.to_csv(
    "processed/kd_pdf_raw_extracted.csv",
    index=False,
    header=False
)

print("\n=== CSV SAVED ===")
print("Total rows:", len(pdf_df))
print("Total columns:", len(pdf_df.columns))
print("Saved: processed/kd_pdf_raw_extracted.csv")