"""Modal wrapper for Tahoe-100M DI computation (script 01).

Run with: modal run --detach experiments/modal_01_tahoe.py

Downloads ~89 GB of parquet shards from HuggingFace, computes DI +
bootstrap CIs for each drug, saves results to Modal volume.

After completion, download results:
    modal volume get di-atlas-results tahoe_di.csv results/tahoe_di.csv
"""
import modal

app = modal.App("di-atlas-tahoe")

vol = modal.Volume.from_name("di-atlas-results", create_if_missing=True)

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install(
        "numpy",
        "scipy",
        "pandas",
        "scikit-learn",
        "matplotlib",
        "pyarrow",
        "huggingface-hub>=1.22.0",
        "tqdm",
    )
)


@app.function(
    image=image,
    volumes={"/results": vol},
    timeout=86400,
    memory=32768,
    cpu=4,
)
def compute_tahoe_di():
    import json
    from datetime import datetime
    from pathlib import Path

    import numpy as np
    import pandas as pd
    from huggingface_hub import HfApi, hf_hub_download
    from scipy import stats
    from sklearn.linear_model import LinearRegression

    DATASET_ID = "tahoebio/Tahoe-100M"
    SUBDIR = "metadata/pseudobulk_differential_expression"
    MIN_CELLS = 20
    MIN_CONTEXTS = 5
    N_BOOTSTRAP = 1000

    ts = lambda: datetime.now().strftime("%H:%M:%S")

    def direction_instability(signatures):
        norms = np.linalg.norm(signatures, axis=1, keepdims=True)
        norms = np.maximum(norms, 1e-10)
        unit = signatures / norms
        cosines = unit @ unit.T
        n = cosines.shape[0]
        triu_idx = np.triu_indices(n, k=1)
        return 1.0 - float(cosines[triu_idx].mean())

    def magnitude_cv(signatures):
        norms = np.linalg.norm(signatures, axis=1)
        if norms.mean() < 1e-10:
            return 0.0
        return float(norms.std() / norms.mean())

    print(f"[{ts()}] Computing Tahoe-100M direction instability")

    api = HfApi()
    files = api.list_repo_tree(DATASET_ID, path_in_repo=SUBDIR, repo_type="dataset")
    shard_paths = [f.rfilename for f in files if f.rfilename.endswith(".parquet")]
    print(f"[{ts()}] Found {len(shard_paths)} parquet shards")

    drug_cellline_genes: dict[tuple[str, str], dict[str, float]] = {}

    for i, shard_path in enumerate(shard_paths):
        if i % 50 == 0:
            print(f"[{ts()}] Processing shard {i+1}/{len(shard_paths)}")

        local = hf_hub_download(
            repo_id=DATASET_ID,
            filename=shard_path,
            repo_type="dataset",
        )
        df = pd.read_parquet(
            local,
            columns=["drug", "Cell_Name_Vevo", "n_cells_trt", "gene_name", "log2FoldChange"],
        )
        df = df[df["n_cells_trt"] >= MIN_CELLS].dropna(subset=["log2FoldChange"])

        for (drug, cell), group in df.groupby(["drug", "Cell_Name_Vevo"]):
            key = (drug, cell)
            if key not in drug_cellline_genes:
                drug_cellline_genes[key] = {}
            drug_cellline_genes[key].update(
                dict(zip(group["gene_name"], group["log2FoldChange"]))
            )

    print(f"[{ts()}] {len(drug_cellline_genes)} drug-cell pairs")

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

    rng = np.random.default_rng(42)
    rows = []
    n_drugs = len(signatures_by_drug)
    for idx, (pid, sigs) in enumerate(signatures_by_drug.items()):
        if idx % 100 == 0:
            print(f"[{ts()}] DI bootstrap: {idx}/{n_drugs}")
            if rows:
                pd.DataFrame(rows).to_csv(Path("/results/tahoe_di_partial.csv"), index=False)
                vol.commit()
        n_ctx = sigs.shape[0]

        di = direction_instability(sigs)
        mcv = magnitude_cv(sigs)
        mean_norm = float(np.linalg.norm(sigs, axis=1).mean())

        boot_dis = []
        for _ in range(N_BOOTSTRAP):
            boot_idx = rng.choice(n_ctx, size=n_ctx, replace=True)
            boot_dis.append(direction_instability(sigs[boot_idx]))
        boot_dis = np.array(boot_dis)

        rows.append({
            "perturbation_id": pid,
            "di_raw": di,
            "di_raw_ci_low": float(np.percentile(boot_dis, 2.5)),
            "di_raw_ci_high": float(np.percentile(boot_dis, 97.5)),
            "magnitude_cv": mcv,
            "mean_norm": mean_norm,
            "n_contexts": n_ctx,
        })

    df = pd.DataFrame(rows)
    print(f"[{ts()}] {len(df)} drugs computed")

    if len(df) > 0:
        X = df["mean_norm"].values.reshape(-1, 1)
        y = df["di_raw"].values
        reg = LinearRegression().fit(X, y)
        pred = reg.predict(X).flatten()
        df["di_corrected"] = (y - pred) + reg.intercept_
        df["di_corrected_ci_low"] = (df["di_raw_ci_low"].values - pred) + reg.intercept_
        df["di_corrected_ci_high"] = (df["di_raw_ci_high"].values - pred) + reg.intercept_

    df = df.rename(columns={"perturbation_id": "drug"})
    df["dataset"] = "Tahoe-100M"

    out_path = Path("/results/tahoe_di.csv")
    df.to_csv(out_path, index=False)
    vol.commit()
    print(f"[{ts()}] Saved {len(df)} drugs to {out_path}")

    return f"Done: {len(df)} drugs saved to Modal volume di-atlas-results/tahoe_di.csv"


@app.local_entrypoint()
def main():
    result = compute_tahoe_di.remote()
    print(result)
