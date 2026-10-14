import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
file_path = "raw/pone.0181585.s004.xlsx"

excel_file = pd.ExcelFile(file_path)

df = pd.read_excel(
    file_path,
    sheet_name = excel_file.sheet_names[0]
)
df_clean = df[
    df["Regno"].notna()
].copy()

              
compounds = df_clean[
    ["Regno", "Compound", "Chemotype", "Smiles"]
].copy()

compounds = compounds.rename(
    columns={
        "Regno": "compound_id", 
        "Compound": "compound_name", 
        "Smiles": "smiles_original"
    }
)
print("\nCompound shape:")
print(compounds.shape)

print("\nCompound columns:")
print(compounds.columns.tolist())

print("\nFirst 10 compounds:")
print(compounds.head(10))

first_smiles = compounds["smiles_original"].iloc[0]

print("\nFirst SMILES:")
print(first_smiles)

first_mol = Chem.MolFromSmiles(first_smiles)

print("\nRDKit molecule:")
print(first_mol)

def parse_smiles(smiles): 
    if pd.isna(smiles):
        return None
    return Chem.MolFromSmiles(smiles)

compounds["mol"] = compounds["smiles_original"].apply(parse_smiles)
print("\nValid molucules:")
print(compounds["mol"].notna().sum())

print("\nInvalid molecules:")
print(compounds["mol"].isna().sum())

compounds["mol"]

compounds["smiles_canonical_rdkit"] = compounds["mol"].apply(
    lambda mol: Chem.MolToSmiles(mol) if mol is not None else None
)

compounds["inchikey_rdkit"] = compounds["mol"].apply(
    lambda mol: Chem.MolToInchiKey(mol) if mol is not None else None
)
print("\nInChKeys:")
print(
    compounds[
        [
        "compound_id", 
        "smiles_canonical_rdkit", 
        "inchikey_rdkit"
    ]
  ].head()  
)

compounds["structure_status"] = compounds["mol"].apply(
    lambda mol: "parsed" if mol is not  None else "failed_to_parse"
)
compounds["notes"] = compounds["mol"].apply(
    lambda mol: "" if mol is not None
    else "Original S4 SMILES could not be parsed by RDKit"
)

print("\nOriginal vs canonical SMILES:")
print(
    compounds[
        ["compound_id","smiles_original", "smiles_canonical_rdkit" ]
    ].head(10)
)
print("\nMissing canonical SMILES:")
print(compounds["smiles_canonical_rdkit"].isna().sum())

print("\nNumber of compounds:")
print(len(compounds))

print("Number of unique canocical smiles_rdkit: ")
print(compounds["smiles_canonical_rdkit"].nunique())
print("Number of duplicated canocical smiles_rdkit: ")
print(compounds["smiles_canonical_rdkit"].duplicated().sum())

duplicate_structures = compounds[
    compounds["smiles_canonical_rdkit"].duplicated(keep=False)
].sort_values("smiles_canonical_rdkit")

print("\nDuplicated strcutures: ")
print(
    duplicate_structures[
        [
        "compound_id", 
        "compound_name", 
        "smiles_original",
        "smiles_canonical_rdkit"
    ]
  ]
)
print(
    duplicate_structures[
        [
        "compound_id", 
        "compound_name", 
        "Chemotype", 
        "smiles_original",
        "smiles_canonical_rdkit"
    ]
   ].to_string(index=False)
)
compounds["molecular_weight"] = compounds["mol"].apply(
    Descriptors.MolWt
)
compounds["logp"] = compounds["mol"].apply(
    Descriptors.MolLogP
)
compounds["h_bond_donors"] = compounds["mol"].apply(
    Descriptors.NumHDonors
)
compounds["h_bond_acceptors"] = compounds["mol"].apply(
    Descriptors.NumHAcceptors
)
print("\nMolecular descriptors:")
print(
    compounds[
        [
        "compound_id", 
        "molecular_weight",
        "logp",
        "h_bond_donors", 
        "h_bond_acceptors"
    ]
  ].head(10) 
)

print("\nMissing descriptor values:")

descriptor_columns = [
    "molecular_weight", 
    "logp", 
    "h_bond_donors", 
    "h_bond_acceptors"
]
print(compounds[descriptor_columns].isna().sum())
print("\nDescriptor summary:")
print(compounds[descriptor_columns].describe())

print("\nMissing important values: ")

important_columns = [
    "compound_id",
    "compound_name",
    "smiles_original",
    "smiles_canonical_rdkit"
]
print(compounds[important_columns].isna().sum())

final_columns = [
    "compound_id",
    "compound_name",
    "Chemotype",
    "smiles_original",
    "smiles_canonical_rdkit",
    "inchikey_rdkit",
    "structure_status",
    "notes",
    "molecular_weight",
    "logp",
    "h_bond_donors",
    "h_bond_acceptors"
]
compounds_final = compounds[final_columns].copy()
print("\nFinal compound shape:")
print(compounds_final.shape)

print("\nFinal columns:")
print(compounds_final.columns.tolist())

print("\nFinal preview: ")
print(compounds_final.head())

compounds_final.to_csv(
    "processed/compounds.csv",
    index=False
)
print("\nSved processed/compounds.csv successfully!")