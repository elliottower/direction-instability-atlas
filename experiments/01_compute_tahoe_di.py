"""Compute direction instability on Tahoe-100M pseudobulk DE signatures.

Reads parquet shards from HuggingFace (1,026 files, ~89 GB total).
Must be run on Modal — too large for local disk.

For each drug: build signature matrix (cell_lines x genes) from
log2FoldChange values, filter cell_lines with n_cells_trt < 20,
compute DI + bootstrap CIs + magnitude correction.

Output: results/tahoe_di.csv
"""
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from huggingface_hub import HfApi, hf_hub_download

from experiments.utils import compute_di_table, magnitude_correct

DATASET_ID = "tahoebio/Tahoe-100M"
SUBDIR = "metadata/pseudobulk_differential_expression"
MIN_CELLS = 20
MIN_CONTEXTS = 5
RESULTS_PATH = Path("results/tahoe_di.csv")


def build_signature_matrices(shard_paths: list[str]) -> dict[str, np.ndarray]:
    """Stream parquet shards and build per-drug signature matrices."""
    drug_cellline_genes: dict[tuple[str, str], dict[str, float]] = {}
    all_genes: set[str] = set()

    ts = lambda: datetime.now().strftime("%H:%M:%S")

    for i, shard_path in enumerate(shard_paths):
        if i % 100 == 0:
            print(f"[{ts()}] Processing shard {i+1}/{len(shard_paths)}")

        local = hf_hub_download(
            repo_id=DATASET_ID,
            filename=shard_path,
            repo_type="dataset",
        )
        df = pd.read_parquet(local)

        for _, row in df.iterrows():
            drug = row["drug"]
            cell = row["Cell_Name_Vevo"]
            n_cells = row.get("n_cells_trt", 0)

            if n_cells < MIN_CELLS:
                continue

            key = (drug, cell)
            if key not in drug_cellline_genes:
                drug_cellline_genes[key] = {}
            drug_cellline_genes[key][row["gene_name"]] = row["log2FoldChange"]
            all_genes.add(row["gene_name"])

    print(f"[{ts()}] {len(drug_cellline_genes)} drug-cell pairs, {len(all_genes)} genes")

    drugs: dict[str, list[tuple[str, dict]]] = {}
    for (drug, cell), genes in drug_cellline_genes.items():
        drugs.setdefault(drug, []).append((cell, genes))

    signatures_by_drug: dict[str, np.ndarray] = {}
    for drug, cell_genes_list in drugs.items():
        if len(cell_genes_list) < MIN_CONTEXTS:
            continue

        gene_sets = [set(genes.keys()) for _, genes in cell_genes_list]
        shared_genes = sorted(set.intersection(*gene_sets))
        if len(shared_genes) == 0:
            continue

        gene_idx = {g: i for i, g in enumerate(shared_genes)}
        sigs = np.zeros((len(cell_genes_list), len(shared_genes)))
        for j, (cell, genes) in enumerate(cell_genes_list):
            for gene in shared_genes:
                sigs[j, gene_idx[gene]] = genes[gene]
        signatures_by_drug[drug] = sigs

    print(f"[{ts()}] {len(signatures_by_drug)} drugs with >= {MIN_CONTEXTS} cell lines")
    return signatures_by_drug


def main():
    ts = lambda: datetime.now().strftime("%H:%M:%S")
    print(f"[{ts()}] Computing Tahoe-100M direction instability")

    api = HfApi()
    files = api.list_repo_tree(DATASET_ID, path_in_repo=SUBDIR, repo_type="dataset")
    shard_paths = [f.rfilename for f in files if f.rfilename.endswith(".parquet")]
    print(f"[{ts()}] Found {len(shard_paths)} parquet shards")

    signatures = build_signature_matrices(shard_paths)

    print(f"[{ts()}] Computing DI + bootstrap CIs (1000 resamples)")
    df = compute_di_table(signatures, n_bootstrap=1000, min_contexts=MIN_CONTEXTS)

    print(f"[{ts()}] Magnitude correction")
    df = magnitude_correct(df)
    df = df.rename(columns={"perturbation_id": "drug"})
    df["dataset"] = "Tahoe-100M"

    RESULTS_PATH.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(RESULTS_PATH, index=False)
    print(f"[{ts()}] Saved {len(df)} drugs to {RESULTS_PATH}")


if __name__ == "__main__":
    main()
