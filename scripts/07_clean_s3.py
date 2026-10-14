
import pandas as pd
import numpy as np
import re
from pathlib import Path

# --------------------------------------------------
# 1. FILE PATHS
# --------------------------------------------------

INPUT = Path("processed/kd_pdf_raw_extracted.csv")
OUTPUT = Path("processed/kd_long_candidate.csv")
REVIEW = Path("processed/kd_needs_review.csv")
SPECIAL = Path("processed/kd_special_construct_review.csv")

df = pd.read_csv(INPUT, header=None, dtype=str).fillna("")

records = []
review = []
special_records = []

# --------------------------------------------------
# 2. HELPER FUNCTIONS
# --------------------------------------------------

def clean(value):
    return re.sub(r"\s+", " ", str(value)).strip()


def parse_kd(text):
    text = clean(text).replace(",", "")

    match = re.fullmatch(
        r"(<=|>=|<|>|=)?\s*(\d+(?:\.\d+)?)",
        text
    )

    if not match:
        return None, None

    operator = match.group(1) or "="
    value = float(match.group(2))

    if value <= 0:
        return None, None

    return operator, value


def parse_target(text):
    text = clean(text)

    match = re.fullmatch(
        r"(.+?)\s+([A-Za-z0-9_-]+)\s*(<=|>=|<|>|=)",
        text
    )

    if match:
        return match.group(1), match.group(2), match.group(3)

    return None, None, None


# --------------------------------------------------
# 3. PROCESS ROWS
# --------------------------------------------------

for index, row in df.iterrows():

    cells = [
        clean(row.iloc[i]) if i < len(row) else ""
        for i in range(7)
    ]

    compound = cells[0]

    # Recover compound numbers split by PDF extraction.
    if compound == "PFE-PKIS":
        match = re.match(r"^(\d+)\b", cells[1])
        if match:
            compound = f"PFE-PKIS {match.group(1)}"

    if not compound or compound.lower() in {
        "compound name", "compound", "nan"
    }:
        continue

    target = None
    gene = None
    operator = None
    kd = None
    special_layout = False

    # Layout B: separate columns.
    if cells[5] in {"=", ">", "<", ">=", "<="}:
        target = cells[3]
        gene = cells[4]
        operator, kd = parse_kd(cells[5] + cells[6])

    # Layout A: combined target and gene.
    else:
        for target_col in [3, 4]:
            if target_col + 1 >= len(cells):
                continue

            t, g, op = parse_target(cells[target_col])

            if op:
                parsed_op, parsed_kd = parse_kd(
                    op + cells[target_col + 1]
                )

                if parsed_kd is not None:
                    target = t
                    gene = g
                    operator = parsed_op
                    kd = parsed_kd
                    break

    # Layout C: a complete kinase construct followed
    # by an operator and Kd value.
    #
    # Preserve the FULL construct as its identifier.
    # Do not automatically map it to a parent kinase.
    if kd is None:
        for col in [3, 4]:
            if col + 1 >= len(cells):
                continue

            match = re.fullmatch(
                r"(.+?)\s*(<=|>=|<|>|=)\s*",
                cells[col]
            )

            if not match:
                continue

            candidate_target = clean(match.group(1))
            candidate_operator = match.group(2)

            parsed_op, parsed_kd = parse_kd(
                candidate_operator + cells[col + 1]
            )

            if (
                parsed_kd is not None
                and re.search(r"[A-Za-z]", candidate_target)
            ):
                target = candidate_target
                gene = candidate_target
                operator = parsed_op
                kd = parsed_kd
                special_layout = True
                break

    if not target or not gene or kd is None:
        review.append({
            "source_row": index,
            "compound_name": compound,
            "raw_columns": repr(cells),
            "reason": "Could not confidently parse row"
        })
        continue

    record = {
        "source_row": index,
        "compound_name": compound,
        "target_original": target,
        "target_id": gene,
        "kd_operator": operator,
        "kd_nm": kd,
        "pkd": 9 - np.log10(kd) if operator == "=" else np.nan,
        "source_file": "pone.0181585.s003.pdf",
        "extraction_method": "pdfplumber"
    }

    records.append(record)

    if special_layout:
        special_records.append({
            **record,
            "raw_columns": repr(cells),
            "verification_status": "needs_source_check"
        })

# --------------------------------------------------
# 4. SAVE OUTPUTS
# --------------------------------------------------

result = pd.DataFrame(records, columns=[
    "source_row",
    "compound_name",
    "target_original",
    "target_id",
    "kd_operator",
    "kd_nm",
    "pkd",
    "source_file",
    "extraction_method"
])

review_df = pd.DataFrame(review, columns=[
    "source_row",
    "compound_name",
    "raw_columns",
    "reason"
])

special_df = pd.DataFrame(special_records, columns=[
    *result.columns,
    "raw_columns",
    "verification_status"
])

OUTPUT.parent.mkdir(parents=True, exist_ok=True)

result.to_csv(OUTPUT, index=False)
review_df.to_csv(REVIEW, index=False)
special_df.to_csv(SPECIAL, index=False)

print("\n=== S3 CLEANING SUMMARY ===")
print("Raw extracted rows:", len(df))
print("Candidate Kd measurements:", len(result))
print("Rows requiring review:", len(review_df))
print("Special-layout records requiring source check:", len(special_df))

print("\n=== KD OPERATORS ===")
print(result["kd_operator"].value_counts())

print("\n=== UNIQUE COMPOUNDS AND TARGETS ===")
print("Compounds:", result["compound_name"].nunique())
print("Targets:", result["target_id"].nunique())

print("\n=== DUPLICATES ===")
print(
    "Duplicate compound-target rows:",
    result.duplicated(
        ["compound_name", "target_id"], keep=False
    ).sum()
)

print("\nSaved:")
print(OUTPUT)
print(REVIEW)
print(SPECIAL)

print("\n=== FIRST 10 SPECIAL RECORDS ===")
if not special_df.empty:
    print(special_df[[
        "source_row",
        "compound_name",
        "target_original",
        "kd_operator",
        "kd_nm"
    ]].head(10).to_string(index=False))
