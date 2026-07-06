"""Compute direction instability on JUMP-CP ORF overexpression profiles.

Same pipeline as CRISPR but different feature space (722 features,
no PCA correction) and different control filtering (pert_type != "trt").

Output: results/orf_di.csv
"""
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from experiments.utils import compute_di_table, magnitude_correct

PROFILES_PATH = Path("data/jump_cp/orf_profiles.parquet")
METADATA_PATH = Path("data/jump_cp/orf.csv.gz")
RESULTS_PATH = Path("results/orf_di.csv")
MIN_CONTEXTS = 5

FEATURE_PREFIX = "X_"


def main():
    ts = lambda: datetime.now().strftime("%H:%M:%S")
    print(f"[{ts()}] Computing JUMP-CP ORF direction instability")

    profiles = pd.read_parquet(PROFILES_PATH)
    metadata = pd.read_csv(METADATA_PATH)
    print(f"[{ts()}] Profiles: {profiles.shape}, Metadata: {metadata.shape}")

    feature_cols = [c for c in profiles.columns if c.startswith(FEATURE_PREFIX)]
    assert feature_cols, f"No feature columns found. Columns: {list(profiles.columns[:20])}"
    print(f"[{ts()}] Using {len(feature_cols)} feature columns")

    jcp_to_gene = dict(zip(metadata["Metadata_JCP2022"], metadata["Metadata_Symbol"]))
    trt_jcps = set(metadata.loc[metadata["Metadata_pert_type"] == "trt", "Metadata_JCP2022"])

    profiles["gene"] = profiles["Metadata_JCP2022"].map(jcp_to_gene)
    profiles = profiles.dropna(subset=["gene"])
    profiles = profiles[profiles["Metadata_JCP2022"].isin(trt_jcps)]
    print(f"[{ts()}] {len(profiles)} wells after filtering to treatment perturbations")

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
    df["dataset"] = "JUMP-CP_ORF"

    RESULTS_PATH.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(RESULTS_PATH, index=False)
    print(f"[{ts()}] Saved {len(df)} genes to {RESULTS_PATH}")


if __name__ == "__main__":
    main()
