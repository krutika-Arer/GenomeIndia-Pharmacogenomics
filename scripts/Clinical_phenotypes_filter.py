#!/usr/bin/env python3

import pandas as pd

# ==========================================================
# Files
# ==========================================================

INPUT_FILE = "Copy of Master_PGx_Final - Master_PGx_Final.tsv"
PHENOTYPE_FILE = "selected_phenotypes.txt"
OUTPUT_FILE = "Master_PGx_SelectedPhenotypes.tsv"

# ==========================================================
# Load phenotype list
# ==========================================================

with open(PHENOTYPE_FILE) as f:
    selected = {
        line.strip().lower()
        for line in f
        if line.strip()
    }

print(f"Loaded {len(selected)} selected phenotypes")

# ==========================================================
# Load master table
# ==========================================================

df = pd.read_csv(
    INPUT_FILE,
    sep="\t",
    dtype=str,
    keep_default_na=False,
    low_memory=False
)

print("Rows before filtering:", len(df))

# ==========================================================
# Function to determine whether to keep a row
# ==========================================================

def keep_row(cell):

    # Keep blank phenotype rows
    if cell is None:
        return True

    cell = str(cell).strip()

    if cell == "":
        return True

    cell_lower = cell.lower()

    # Exact match
    if cell_lower in selected:
        return True

    # Split into individual phenotypes
    phenotypes = [x.strip().lower() for x in cell.split(",")]

    # Keep if ANY phenotype matches
    for p in phenotypes:
        if p in selected:
            return True

    return False

# ==========================================================
# Filter
# ==========================================================

filtered = df[df["Clin_phenotypes"].apply(keep_row)]

print("Rows after filtering:", len(filtered))

# ==========================================================
# Save
# ==========================================================

filtered.to_csv(
    OUTPUT_FILE,
    sep="\t",
    index=False
)

print("\nDone!")
print("Saved as:", OUTPUT_FILE)
