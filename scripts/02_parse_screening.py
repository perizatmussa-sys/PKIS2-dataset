import pandas as pd

file_path = "raw/pone.0181585.s004.xlsx"

excel_file = pd.ExcelFile(file_path)

df = pd.read_excel(
    file_path,
    sheet_name = excel_file.sheet_names[0]
)
print("Raw shape:", df.shape)
print("Rows with Regno:", df["Regno"].notna().sum())
print("Rows without Regno:", df["Regno"].isna().sum())
df_clean = df[df["Regno"].notna()].copy()
print("Clean shape:", df_clean.shape)
#columns that are not individula kinase assay measuremnts 
compound_metadata_columns = [
    "Regno", 
    "Compound", 
    "Chemotype", 
    "Smiles"
]

summary_columns = [
    ">90", 
    ">80", 
    ">70"
]
non_asssay_columns = compound_metadata_columns + summary_columns 
assay_columns = [
    col for col in df_clean.columns
    if col not in non_asssay_columns
]

print("Number of non-assay columns: ", len(non_asssay_columns))
print("Number of potential assay columns: ", len(assay_columns))

print("\nFirst 10 potential assay columns: ")
print(assay_columns[:10])

print("\nLast 10 potential assay columns: ")
print(assay_columns[-10:])

screening_long = df_clean.melt(
    id_vars=["Regno", "Compound"], 
    value_vars = assay_columns, 
    var_name = "target_original", 
    value_name="raw_value"
)

screening_long = screening_long.rename(
    columns={
        "Regno": "compound_id",
        "Compound": "compound_name"
    }
)
screening_long["raw_measurement_type"] = "percent_inhibition"
screening_long["standardized_value"] = screening_long["raw_value"]
screening_long["standardized_measurement_type"] = "percent_inhibition"

pd.set_option("display.max_columns", None)
pd.set_option("display.width", None)

print("\nUpdated columns:")
print(screening_long.columns.tolist())

print("\nUpdated first 5 rows: ")
print(screening_long.head())

screening_long["source_file"] = "pone.0181585.s004.xlsx"
screening_long["source_sheet"] = excel_file.sheet_names[0]

print("\nColumns after adding source information:")
print(screening_long.columns.tolist())

print("\nSource information:")
print(
    screening_long[
        ["source_file", "source_sheet"]
    ].head()
)
screening_long["screening_concentration_uM"] = 1.0
print("\nScreening concentration:")
print(screening_long["screening_concentration_uM"].value_counts())

targets = pd.DataFrame({
    "target_original": screening_long["target_original"].unique()
})

print("\nNumber of unique targets:")
print(len(targets))

print("\nFirst 20 targets:")
print(targets.head(20))

targets["target_id"] = targets["target_original"]
targets["species"] = pd.NA
targets["mutation"] = pd.NA
targets["phosphorylation_state"] = pd.NA

targets["target_class"] = "unresolved"

targets["notes"] = (
    "Target metadata do not yet resolved from source documentation"
)

print("Target metadata table:")
print(targets.head(30))
print("\nNumber of targets:")
print(len(targets))

print("\nTarget class counts:")
print(targets["target_class"].value_counts())

screening_long = screening_long.merge(
    targets[["target_original", "target_id"]],
    on = "target_original",
    how = "left"
)

special_keywords = [
    "phosphorylated", 
    "nonphosphorylated", 
    "mutant",
    "human",
    "mouse",
    "rat"
]
for keyword in special_keywords:
    matches = targets[
        targets["target_original"].str.contains(
            keyword, 
            case=False, 
            na=False
        )
    ]
    print(f"\nTargets containing '{keyword}':")
    print(matches["target_original"].tolist())

print("\nTargets with IDs:")
print(targets.head(20))

print("\nMissing target IDs in screening data:")
print(screening_long["target_id"].isna().sum())

print("\nScreening shape after merge:")
print(screening_long.shape)

print("\nAll target names:")
for target in targets["target_original"]:
    print(target)

targets["target_class"] = "unresolved"

print("\nTargets before classification:")
print(targets.head(20))

print("\nTarget class counts:")
print(targets["target_class"].value_counts())

screening_long = screening_long.merge(
    targets[["target_id", "target_class"]], 
    on = "target_id",
    how="left"
)
print("\nMissing target classes:")
print(screening_long["target_class"].isna().sum())

print("\nTarget class counts:")
print(screening_long["target_class"].value_counts())

print("\nShape after adding target class:")
print(screening_long.shape)

print("\nFinal screening columns:")
for col in screening_long.columns:
    print(col)

final_columns = [
    "compound_id",
    "compound_name",
    "target_id",
    "target_original",
    "target_class",
    "raw_value",
    "raw_measurement_type",
    "standardized_value",
    "standardized_measurement_type",
    "screening_concentration_uM",
    "source_file",
    "source_sheet"
]

screening_long = screening_long[final_columns]

print("\nFinal shape:")
print(screening_long.shape)

print("\nFinal columns:")
print(screening_long.columns.tolist())

print("\nFinal preview:")
print(screening_long.head())
targets["phosphorylation_state"] = "not_specified"

# -----------------------------------------
# Record special assay annotations
# -----------------------------------------

# Default: no special annotation detected
targets["notes"] = "none"

# Autoinhibited assays
mask = targets["target_original"].str.contains(
    "autoinhibited",
    case=False,
    na=False
)
targets.loc[mask, "notes"] = "autoinhibited assay"

# Cyclin-associated assays
mask = targets["target_original"].str.contains(
    "cyclin",
    case=False,
    na=False
)
targets.loc[mask, "notes"] = "cyclin-associated assay"

# JH1 catalytic domain assays
mask = targets["target_original"].str.contains(
    "JH1domain-catalytic",
    case=False,
    na=False
)
targets.loc[mask, "notes"] = "JH1 catalytic domain assay"

# JH2 pseudokinase domain assays
mask = targets["target_original"].str.contains(
    "JH2domain-pseudokinase",
    case=False,
    na=False
)
targets.loc[mask, "notes"] = "JH2 pseudokinase domain assay"

# Kinase domain 1
mask = targets["target_original"].str.contains(
    "Kin.Dom.1",
    case=False,
    regex=False,
    na=False
)
targets.loc[mask, "notes"] = "kinase domain 1 assay"

# Kinase domain 2
mask = targets["target_original"].str.contains(
    "Kin.Dom.2",
    case=False,
    regex=False,
    na=False
)
targets.loc[mask, "notes"] = "kinase domain 2 assay"

print("\nAssay annotation counts:")
print(targets["notes"].value_counts(dropna=False))

print("\nTargets with special assay annotations:")
print(
    targets.loc[
        targets["notes"] != "none",
        ["target_original", "notes"]
    ].to_string(index=False)
)

targets.loc[
    targets["target_original"].str.contains(
        "nonphosphorylated",
        case=False,
        na=False
    ),
    "phosphorylation_state"
] = "nonphosphorylated"

targets.loc[
    targets["target_original"].str.contains(
        r"(?<!non)phosphorylated",
        case=False,
        regex=True,
        na=False
    ),
    "phosphorylation_state"
] = "phosphorylated"

print("\nPhosphorylation-state classification:")
print(
    targets["phosphorylation_state"]
    .value_counts(dropna=False)
)
mutation_pattern = r"(?<![A-Za-z0-9])[A-Z]\d+[A-Z](?![A-Za-z0-9])"

mutation_mask = targets["target_original"].str.contains(
    mutation_pattern,
    regex=True,
    na=False
)

targets.loc[
    mutation_mask,
    "mutation"
] = "mutant"

print("\nMutation classification:")
print(targets["mutation"].value_counts(dropna=False))

print("\nTargets classified as mutant:")
print(
    targets.loc[
        targets["mutation"] == "mutant",
        ["target_original", "target_id"]
    ].to_string(index=False)
)
# Save the final screening table
screening_long.to_csv(
    "processed/screening_long.csv",
    index=False
)

# Save the target table
targets.to_csv(
    "processed/targets.csv",
    index=False
)

print("\nFiles saved successfully!\n")

print("Targets containing parentheses:")
parenthesis_targets = targets[
    targets["target_original"].str.contains(
        r"\(",
        regex=True,
        na=False
    )
]
print(
    parenthesis_targets[
        ["target_original", "target_id"]
    ].to_string(index=False)
)
print("\nNumber of targets containing parentheses:")
print(len(parenthesis_targets))

print("\nAll speacially annotated targets:")

for number, target in enumerate(
    parenthesis_targets["target_original"],
    start=1
):
    print(number, target)

print("\nTarget class values:")
print(targets["target_class"].value_counts(dropna=False))

print("\nNotes values:")
print(targets["notes"].value_counts(dropna=False).head(20))

species_terms = [
    "mouse", "rat", "bovine", "canine",
    "porcine","zebrafish", "chimp",
    "human"
]

species_like = targets[
    targets["target_original"]
    .str.lower()
    .str.contains("|".join(species_terms),na=False)
]

print("\nTargets containing obvious species terms:")
print(
    species_like[
        ["target_original", "target_id"]
    ].to_string(index=False)
)

print("\nNumber with obvious species terms:")
print(len(species_like))
special_mask = (
    targets["target_original"].str.contains(
        r"\(|phosphorylated|nonphosphorylated|autoinhibited|cyclin",
        case=False,
        regex=True,
        na=False
    )
)

special_targets = targets.loc[
    special_mask,
    ["target_original", "target_id"]
]

print("\nSpecially annotated assay targets:")
print(special_targets.to_string(index=False))

print("\nNumber of specially annotated assay targets:")
print(len(special_targets))

# -----------------------------------------
# Inspect candidate variant assays
# -----------------------------------------

targets["mutation"] = "not_specified"

# Detect explicit amino-acid substitutions such as S808G
mutation_pattern = r"(?<![A-Za-z0-9])[A-Z]\d+[A-Z](?![A-Za-z0-9])"

mutation_mask = targets["target_original"].str.contains(
    mutation_pattern,
    regex=True,
    na=False
)

targets.loc[mutation_mask, "mutation"] = "mutant"

print("\nMutation classification:")
print(targets["mutation"].value_counts(dropna=False))

print("\nTargets classified as mutant:")
print(
    targets.loc[
        targets["mutation"] == "mutant",
        ["target_original", "target_id"]
    ].to_string(index=False)
)


# -----------------------------------------
# 2. THEN identify candidate special assays
# -----------------------------------------

candidate_variant_mask = (
    targets["target_original"].str.contains(
        r"phosphorylated|nonphosphorylated|autoinhibited|"
        r"cyclin|Kin\.Dom|JH1domain|JH2domain",
        case=False,
        regex=True,
        na=False
    )
)

candidate_variants = targets.loc[
    candidate_variant_mask,
    [
        "target_original",
        "target_id",
        "mutation",
        "phosphorylation_state",
        "notes"
    ]
].copy()


# -----------------------------------------
# 3. Print candidate assays
# -----------------------------------------

print("\nCandidate variant/special-form assays:")
print(candidate_variants.to_string(index=False))

print("\nNumber of candidate variant/special-form assays:")
print(len(candidate_variants))


# -----------------------------------------
# 4. Breakdown
# -----------------------------------------

print("\n--- BREAKDOWN OF CANDIDATE ASSAYS ---")

for _, row in candidate_variants.iterrows():
    print(
        row["target_original"],
        "| mutation:", row["mutation"],
        "| phospho:", row["phosphorylation_state"],
        "| notes:", row["notes"]
    )


# -----------------------------------------
# 5. Summary counts
# -----------------------------------------

print("\nCandidate counts by annotation:")
print(candidate_variants["notes"].value_counts(dropna=False))

print("\nCandidate counts by mutation:")
print(candidate_variants["mutation"].value_counts(dropna=False))

print("\nCandidate counts by phosphorylation state:")
print(
    candidate_variants["phosphorylation_state"]
    .value_counts(dropna=False)
)


# -----------------------------------------
# 6. Final candidate table
# -----------------------------------------

print("\n=== ALL CANDIDATES ===")

print(
    candidate_variants[
        [
            "target_original",
            "mutation",
            "phosphorylation_state",
            "notes"
        ]
    ].to_string(index=False)
)
import re 

def get_parent_target(name):
    name = str(name)

    parent = re.sub(r"\(.*?\)", "", name)

    parent = re.sub(
        r"-(nonphosphorylated|phosphorylated|autoinhibited|cyclinD1|cyclinD3)$",
        "",
        parent,
        flags=re.IGNORECASE 
    )

    return parent.strip()

targets["parent_target"] = (
    targets["target_original"]
    .apply(get_parent_target)
)

parent_counts = targets["parent_target"].value_counts()

multiple_forms = parent_counts[parent_counts > 1]

print("\nParents having multiple assay forms:")
print(multiple_forms)

print("\nAssays belonging to those parents:")

multi_assays = targets[
    targets["parent_target"].isin(multiple_forms.index)
][
     ["parent_target", "target_original",
     "mutation", "phosphorylation_state", "notes"]
].sort_values(["parent_target", "target_original"])

print(multi_assays.to_string(index=False))

print("\nNumber of assay rows involved:")
print(len(multi_assays))

print("\n=== MULTIPLE-FORM ASSAYS ONLY ===")

print(
    multi_assays[
        [
            "parent_target",
            "target_original",
            "mutation",
            "phosphorylation_state",
            "notes"
        ]
    ].to_string(index=False)
)

variant_targets = [
    
]
targets["species"] = "human"

# Start by assigning WT
targets["wild_type_or_mutant"] = "wild_type"

# Change the known variants
targets.loc[
    targets["target_original"].isin(variant_targets),
    "wild_type_or_mutant"
] = "variant"

# Check result
print("\nWT / variant classification:")
print(targets["wild_type_or_mutant"].value_counts())

# -----------------------------------------
# Check unique biological kinase identities
# -----------------------------------------

print("\n=== TARGET IDENTITY CHECK ===")

print("Total S4 assay columns:")
print(len(targets))

print("\nUnique parent kinase identities:")
print(targets["parent_target"].nunique())

print("\nDifference between assay columns and parent identities:")
print(len(targets) - targets["parent_target"].nunique())

# -----------------------------------------
# Final target metadata
# -----------------------------------------

# The S4 panel is the human kinase panel
targets["species"] = "human"

# Biological parent kinase
targets["target_name_clean"] = targets["parent_target"]

# Explicit mutation status
targets["wild_type_or_mutant"] = "wild_type"

targets.loc[
    targets["mutation"] == "mutant",
    "wild_type_or_mutant"
] = "mutant"

# Parent target for assay forms
targets["parent_target_if_mutant"] = pd.NA

targets.loc[
    targets["wild_type_or_mutant"] == "mutant",
    "parent_target_if_mutant"
] = targets["parent_target"]

# -----------------------------------------
# SAVE FINAL TARGET METADATA
# -----------------------------------------

targets.to_csv(
    "processed/targets.csv",
    index=False
)

print("\nSaved final target metadata:")
print("processed/targets.csv")

print("\nFinal saved columns:")
print(targets.columns.tolist())

print("\n=== FINAL TARGET METADATA CHECK ===")

print("Assay rows:", len(targets))
print("Unique kinase identities:", targets["target_name_clean"].nunique())

print("\nSpecies:")
print(targets["species"].value_counts(dropna=False))

print("\nWT / mutant:")
print(targets["wild_type_or_mutant"].value_counts(dropna=False))