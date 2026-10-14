
import pandas as pd
from pathlib import Path

folder = Path("processed")

s3 = pd.read_csv(folder / "kd_long_candidate.csv", dtype=str)
s4 = pd.read_csv(folder / "targets.csv", dtype=str)
unresolved = pd.read_csv(
    folder / "s3_unresolved_targets.csv",
    dtype=str
)

def normalize(value):
    if pd.isna(value):
        return ""
    return str(value).strip().upper()

# Available S4 assay names
s4_names = {}

for _, row in s4.iterrows():
    assay_id = row["target_id"]

    for column in ["target_id", "target_original",
                   "target_name_clean"]:
        name = normalize(row[column])

        if name:
            s4_names.setdefault(name, set()).add(assay_id)

results = []

for gene in unresolved["target_id"]:
    gene = normalize(gene)

    # Find all original kinase names associated
    # with this gene symbol in S3.
    original_names = set(
        s3.loc[
            s3["target_id"].map(normalize) == gene,
            "target_original"
        ].dropna().map(normalize)
    )

    matches = set()

    for name in original_names:
        matches.update(s4_names.get(name, set()))

    if len(matches) == 1:
        status = "candidate_unique"
    elif len(matches) > 1:
        status = "ambiguous_multiple"
    else:
        status = "unresolved"

    results.append({
        "s3_gene": gene,
        "s3_original_names": "; ".join(sorted(original_names)),
        "s4_candidates": "; ".join(sorted(matches)),
        "status": status
    })

result = pd.DataFrame(results)

result.to_csv(
    folder / "s3_target_mapping_candidates.csv",
    index=False
)

print("\n=== UPDATED TARGET MATCHING ===")
print(result["status"].value_counts())

print("\n=== FIRST 20 RESULTS ===")
print(result.head(20).to_string(index=False))
