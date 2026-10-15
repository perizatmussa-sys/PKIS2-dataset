
from pathlib import Path
import pandas as pd
import numpy as np

# Locate the project
root = Path(__file__).resolve().parents[1]

input_file = root / "processed" / "screening_long.csv"
output_dir = root / "splits_25"
output_dir.mkdir(exist_ok=True)

# Load the PKIS2 dataset
df = pd.read_csv(input_file)

# Check the expected number of observations
assert len(df) == 261870, "Unexpected dataset size!"

# Randomly shuffle rows, reproducibly
rng = np.random.default_rng(42)
order = rng.permutation(len(df))
df = df.iloc[order].reset_index(drop=True)

# Create 25 separate training datasets
for i in range(25):
    start = i * 10000
    end = start + 10000

    train = df.iloc[start:end]

    filename = f"train_{i+1:02d}.csv"
    train.to_csv(output_dir / filename, index=False)

    print(f"{filename}: {len(train)} rows")

# Save all remaining observations for testing
test = df.iloc[250000:]
test.to_csv(output_dir / "test.csv", index=False)

print("\nTest rows:", len(test))
print("Total observations:", len(df))
print("Saved to:", output_dir)
