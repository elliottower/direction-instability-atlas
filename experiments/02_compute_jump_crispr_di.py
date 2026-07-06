"""Compute direction instability on JUMP-CP CRISPR knockout profiles.

Loads well-level profiles (post-sphering, harmony, PCA corrected),
maps JCP2022 IDs to gene symbols, excludes controls, groups by gene,
computes DI + bootstrap CIs + magnitude correction.

Output: results/crispr_di.csv
"""
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from experiments.utils import compute_di_table, magnitude_correct

PROFILES_PATH = Path("data/jump_cp/crispr_profiles.parquet")
METADATA_PATH = Path("data/jump_cp/crispr.csv.gz")
RESULTS_PATH = Path("results/crispr_di.csv")
MIN_CONTEXTS = 5

FEATURE_PREFIX = "X_"
CONTROL_GENES = {"", "non-targeting", "no-guide"}


def main():
    ts = lambda: datetime.now().strftime("%H:%M:%S")
    print(f"[{ts()}] Computing JUMP-CP CRISPR direction instability")

    profiles = pd.read_parquet(PROFILES_PATH)
    metadata = pd.read_csv(METADATA_PATH)
    print(f"[{ts()}] Profiles: {profiles.shape}, Metadata: {metadata.shape}")

    feature_cols = [c for c in profiles.columns if c.startswith(FEATURE_PREFIX)]
    if not feature_cols:
        feature_cols = [c for c in profiles.columns if c.startswith("PC_")]
    assert feature_cols, f"No feature columns found. Columns: {list(profiles.columns[:20])}"
    print(f"[{ts()}] Using {len(feature_cols)} feature columns")

    jcp_to_gene = dict(zip(metadata["Metadata_JCP2022"], metadata["Metadata_Symbol"]))

    profiles["gene"] = profiles["Metadata_JCP2022"].map(jcp_to_gene)
    profiles = profiles.dropna(subset=["gene"])
    profiles = profiles[~profiles["gene"].isin(CONTROL_GENES)]
    print(f"[{ts()}] {len(profiles)} wells after removing controls")

    source_nunique = profiles["Metadata_Source"].nunique() if "Metadata_Source" in profiles.columns else 0
    plate_col = "Metadata_Source" if source_nunique > 1 else "Metadata_Plate"
    signatures_by_gene: dict[str, np.ndarray] = {}
    for gene, group in profiles.groupby("gene"):
        plate_means = group.groupby(plate_col)[feature_cols].mean()
        sigs = plate_means.values.astype(np.float64)
        if sigs.shape[0] >= MIN_CONTEXTS:
            signatures_by_gene[gene] = sigs

    print(f"[{ts()}] {len(signatures_by_gene)} genes with >= {MIN_CONTEXTS} replicates")

    print(f"[{ts()}] Computing DI + bootstrap CIs")
    df = compute_di_table(signatures_by_gene, n_bootstrap=1000, min_contexts=MIN_CONTEXTS)

    print(f"[{ts()}] Magnitude correction")
    df = magnitude_correct(df)
    df = df.rename(columns={"perturbation_id": "gene"})
    df["dataset"] = "JUMP-CP_CRISPR"

    RESULTS_PATH.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(RESULTS_PATH, index=False)
    print(f"[{ts()}] Saved {len(df)} genes to {RESULTS_PATH}")


if __name__ == "__main__":
    main()
