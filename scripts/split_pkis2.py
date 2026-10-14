
from pathlib import Path
import pandas as pd
import numpy as np

# Project folders
root = Path(__file__).resolve().parents[1]
input_file = root / "processed" / "screening_long.csv"
output_dir = root / "splits"
output_dir.mkdir(exist_ok=True)

# Load data
df = pd.read_csv(input_file)

# Shuffle unique compounds reproducibly
compounds = df["compound_id"].dropna().unique()
rng = np.random.default_rng(42)
rng.shuffle(compounds)

# Divide compounds: 80% / 10% / 10%
n_train = int(len(compounds) * 0.8)
n_val = int(len(compounds) * 0.1)

train_ids = set(compounds[:n_train])
val_ids = set(compounds[n_train:n_train + n_val])
test_ids = set(compounds[n_train + n_val:])

train = df[df["compound_id"].isin(train_ids)]
val = df[df["compound_id"].isin(val_ids)]
test = df[df["compound_id"].isin(test_ids)]

# Save CSV files
train.to_csv(output_dir / "train.csv", index=False)
val.to_csv(output_dir / "validation.csv", index=False)
test.to_csv(output_dir / "test.csv", index=False)

# Check no compound appears in multiple splits
assert train_ids.isdisjoint(val_ids)
assert train_ids.isdisjoint(test_ids)
assert val_ids.isdisjoint(test_ids)

print("Training rows:", len(train))
print("Validation rows:", len(val))
print("Testing rows:", len(test))
print("Total:", len(train) + len(val) + len(test))
print("Files saved in:", output_dir)
