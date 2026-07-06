"""Compute viability coefficient of variation from PRISM 24Q2.

PRISM viability is 1D per cell line — DI (cosine) is undefined.
Instead compute CV of LFC across cell lines for each drug.

Input: data/prism/Repurposing_Public_24Q2_LFC_COLLAPSED.csv
Output: results/prism_cv.csv
"""
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

LFC_PATH = Path("data/prism/Repurposing_Public_24Q2_LFC_COLLAPSED.csv")
COMPOUND_PATH = Path("data/prism/Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv")
RESULTS_PATH = Path("results/prism_cv.csv")
MIN_CELL_LINES = 5


def main():
    ts = lambda: datetime.now().strftime("%H:%M:%S")
    print(f"[{ts()}] Computing PRISM viability CV")

    lfc = pd.read_csv(LFC_PATH)
    lfc = lfc[~lfc["broad_id"].str.contains("QC Failure", na=False)]
    print(f"[{ts()}] LFC data: {lfc.shape}")

    compound_meta = pd.read_csv(COMPOUND_PATH)
    compound_meta["broad_id"] = compound_meta["IDs"].str.replace("BRD:", "", regex=False)
    brd_to_name = dict(zip(compound_meta["broad_id"], compound_meta["Drug.Name"]))

    drug_stats = []
    for broad_id, group in lfc.groupby("broad_id"):
        mean_lfc_per_cell = group.groupby("row_id")["LFC"].mean().dropna()
        lfc_values = mean_lfc_per_cell.values
        n_cells = len(lfc_values)

        if n_cells < MIN_CELL_LINES:
            continue

        mean_val = np.mean(lfc_values)
        std_val = np.std(lfc_values, ddof=1)

        if abs(mean_val) < 0.01:
            continue
        cv = std_val / abs(mean_val)

        drug_name = brd_to_name.get(broad_id, broad_id)

        drug_stats.append({
            "broad_id": broad_id,
            "drug_name": drug_name,
            "viability_cv": cv,
            "mean_lfc": mean_val,
            "std_lfc": std_val,
            "n_cell_lines": n_cells,
        })

    df = pd.DataFrame(drug_stats)
    print(f"[{ts()}] {len(df)} drugs with >= {MIN_CELL_LINES} cell lines")

    RESULTS_PATH.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(RESULTS_PATH, index=False)
    print(f"[{ts()}] Saved to {RESULTS_PATH}")


if __name__ == "__main__":
    main()
