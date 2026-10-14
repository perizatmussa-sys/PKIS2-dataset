import pandas as pd
screening_path = "processed/screening_long.csv"
kd_path = "processed/kd_long.csv"

screening_df = pd.read_csv(screening_path)
kd_df = pd.read_csv(kd_path)

print("\nScreening shape:")
print(screening_df.shape)

print("\nKd shape:")
print(kd_df.shape)

print("\nScreening columns:")
print(screening_df.columns.tolist())

print("\nKd columns:")
print(kd_df.columns.tolist())

print("\nUnique compounds in screening:")
print(screening_df["compound_name"].nunique())

print("\nUnique compounds in Kd:")
print(kd_df["Compound Name"].nunique())

screening_compounds = set(screening_df["compound_name"].dropna())
kd_compounds = set(kd_df["Compound Name"].dropna())
common_compounds = screening_compounds & kd_compounds
print("\nCompounds found in both datasets:")
print(len(common_compounds))

print("\nKd compounds not found in screening:")
print(len(kd_compounds - screening_compounds))

unmatched_kd = kd_compounds - screening_compounds

print("\nKd compounds not matched to screening:")

for compound in sorted(unmatched_kd):
    print(compound)
for compound in sorted(unmatched_kd):
    print("\nLooking for:", compound)

    matches = screening_df[
        screening_df["compound_name"]
        .astype(str)
        .str.contains(str(compound),case=False, regex =False)
    ]["compound_name"].unique()

    print("Possible matches:", matches[:10])

kd_df = kd_df.rename(columns={
    "Compound Name": "compound_name",
    "Entrez Gene\nSymbol": "target_id",
    "Kd (nM)" : "Kd_nm"
})
print("\nKd columns after renaming:")
print(kd_df.columns.tolist())

screening_pairs = set(
    zip(
        screening_df["compound_name"],
        screening_df["target_id"]
    )
)

kd_pairs = set(
    zip(
        kd_df["compound_name"],
        kd_df["target_id"]
    )
)
common_pairs = screening_pairs & kd_pairs

print("\nUnique compound-target pairs in screening:")
print(len(screening_pairs))

print("\nUnique compound-target pairs in Kd:")
print(len(kd_pairs))

print("\nCompound-target pairs found in BOTH:")
print(len(common_pairs))

linked_df = screening_df.merge(
    kd_df,
    on=["compound_name", "target_id"],
    how="inner"
)

print("\nLinked dataset shape:")
print(linked_df.shape)
print("\nLnked columns:")
print(linked_df.columns.tolist())
print("\nFirst 10 linked rows:")
print(linked_df.head(10))

print("\FINAL MERGE CHECK")
print("\nLinked dataset shape:")
print(linked_df.shape)

print("\nUnique compound-target pairs in linked data:")
print(
    linked_df[["compound_name","target_id"]]
    .drop_duplicates()
    .shape[0]
)

print("\nMissing Kd values:")
print(linked_df["Kd_nm"].isna().sum())

print("\nMissing compound names:")
print(linked_df["compound_name"].isna().sum())

print("\nMissing target IDs:")
print(linked_df["target_id"].isna().sum())

duplicates = linked_df[
    linked_df.duplicated(
        subset=["compound_name", "target_id"],
        keep=False
    )
]

print("\nRepeated compound-target pairs:")
print(duplicates[
    ["compound_name", "target_id","Kd_nm"]
])

print("\nFull repeated rows:")
print(duplicates.to_string)

output_path = "processed/linked_creening_kd.csv"

linked_df.to_csv(
    output_path, 
    index=False
)

print("\nSaved linked dataset to:")
print(output_path)

screening_targets = set(screening_df["target_id"].dropna())
kd_targets = set(kd_df["target_id"].dropna())
common_targets = screening_targets & kd_targets

print("\nUnique S4 targets:")
print(len(screening_targets))

print("\nUnique S3 targets:")
print(len(kd_targets))

print("\nS3 targets found in S4:")
print(len(common_targets))

print("\nS3 targets not found in S4:")
print(len(kd_targets - screening_targets))

unmatched_targets = kd_targets - screening_targets

print("\nS3 targets not matched to S4:")

for target in sorted(unmatched_targets):
    print(target)

print("\nPOSSIBLE MATCHES FOR UNMATCHED TARGETS:")

for targets in sorted(unmatched_targets):
    possible_matches = [
        s4_target
        for s4_target in screening_targets
        if str(target).lower() in str(s4_target).lower()
        or str(s4_target).lower() in str(target).lower()
    ]

print("\nS3 target:",target)
print("Possible S4 matches:", possible_matches)
possible_target_matches = {}
still_unresolved_targets = []

for target in sorted(unmatched_targets):
    matches = [
        s4_target
        for s4_target in screening_targets
        if str(target).lower() in str(s4_target).lower()
        or str(s4_target).lower() in str(target).lower()
    ]

    if matches:
        possible_target_matches[target] = matches
    else:
        still_unresolved_targets.append(target)

print("\nUnmatched S3 targets:")
print(len(unmatched_targets))

print("\nTargets with possible S4 name macthes:")
print(len(possible_target_matches))

print("\nTargets with NO possible S4 match:")
print(len(still_unresolved_targets))

print("\nPossible target mappings:")

for target, matches in possible_target_matches.items():
    print(target, "->", matches)

mapping_report = {
    "s3_unique_compounds": len(kd_compounds),
    "s3_compounds_matched_to_s4":len(common_compounds),
    "s3_compounds_unresolved":len(unmatched_kd),

    "s3_unique_targets" : len(kd_targets),
    "s3_targets_exactly_matched_to_s4" : len(common_targets),
    "s3_targets_unresolved_exact_match" : len(unmatched_targets),

    "common_compound_target_pairs":len(common_pairs)

}

print("\nSection 8 mapping report")
for key, value in mapping_report.items():
    print(key, ":", value)

unresolved_compounds_df = pd.DataFrame({
    "unresolved_compound": sorted(unmatched_kd)
})

unresolved_compounds_df.to_csv(
    "processed/unresolved_compounds.csv",
    index=False
)

unresolved_targets_df = pd.DataFrame({
    "unresolved_compound": sorted(unmatched_targets)
})

unresolved_targets_df.to_csv(
    "processed/unresolved_targets.csv",
    index=False
)
print("\nSaved unresolved compounds:")
print(len(unresolved_compounds_df))

print("\nSaved unresolved targets:")
print(len(unresolved_targets_df))

duplicate_kd_mappings = kd_df[
    kd_df.duplicated(
        subset=["compound_name","target_id"],
        keep=False
    )
].sort_values(
    ["compound_name","target_id"]
)

print("\nDuplicate S3 compound-target mappings:")

print(duplicate_kd_mappings[
    ["compound_name", "target_id","Kd_nm"]
])

print("\nDuplicate S3 compound-target mappings:")
print(duplicate_kd_mappings[
    ["compound_name", "target_id", "Kd_nm"]
])

print("\nNumber of rows involved in duplciate mappings:")
print(len(duplicate_kd_mappings))

duplicate_kd_mappings.to_csv(
    "processed/duplciate_kd_mappings.csv",
    index=False
)

print("\nDuplciate mappings saved:")
print(len(duplicate_kd_mappings))