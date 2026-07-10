"""Experiment G: cpg0000 pilot cross-cell-line CRISPR essentiality.

Tests whether the essentiality sign-flip (Batch 1 Exp C: rho=-0.33)
replicates in direction when DI is computed across cell lines (A549
and U2OS) rather than plate replicates.

cpg0000 is the JUMP pilot dataset with CRISPR profiles in two cell
lines. With K=2, DI degenerates to 1 - cos(s_A549, s_U2OS). This is
a directional convergent check, not a confirmatory test.

Pre-registered in PREREGISTRATION_BATCH3.md (Batch 3, reviewer-prompted).

Inputs:
  - cpg0000 pilot: Cell Painting Gallery (AWS S3, source_4)
  - Cell line metadata: jump-cellpainting/2024_Chandrasekaran_NatureMethods
  - DepMap Chronos 24Q4: data/depmap/24q4/CRISPRGeneEffect.csv

Output: results/cpg0000_pilot.json
"""
import gzip
import urllib.request
from datetime import datetime
from io import BytesIO
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from experiments.utils import save_results

DEPMAP_PATH = Path("data/depmap/24q4/CRISPRGeneEffect.csv")
CPG0000_CACHE = Path("data/cpg0000_pilot")
RESULTS_PATH = Path("results/cpg0000_pilot.json")

CHRONOS_THRESHOLD = -0.5
ALPHA = 0.0167
CELL_LINES = ["A549", "U2OS"]

S3_BASE = "https://cellpainting-gallery.s3.amazonaws.com/cpg0000-jump-pilot/source_4/workspace/profiles/2020_11_04_CPJUMP1"
METADATA_URL = "https://raw.githubusercontent.com/jump-cellpainting/2024_Chandrasekaran_NatureMethods/main/benchmark/output/experiment-metadata.tsv"
PROFILE_SUFFIX = "_normalized_feature_select_negcon_batch.csv.gz"


def get_crispr_plate_cell_lines():
    """Get plate-to-cell-line mapping from companion metadata."""
    ts = lambda: datetime.now().strftime("%H:%M:%S")
    meta_path = CPG0000_CACHE / "experiment_metadata.tsv"

    if not meta_path.exists():
        print(f"[{ts()}] Downloading experiment metadata...")
        CPG0000_CACHE.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(METADATA_URL, meta_path)

    meta = pd.read_csv(meta_path, sep="\t")
    crispr = meta[
        (meta["Perturbation"] == "crispr")
        & (meta["Batch"] == "2020_11_04_CPJUMP1")
    ]

    plate_to_cellline = dict(
        zip(crispr["Assay_Plate_Barcode"], crispr["Cell_type"])
    )
    print(f"[{ts()}] CRISPR plates: {len(plate_to_cellline)}")
    for cl in CELL_LINES:
        n = sum(1 for v in plate_to_cellline.values() if v == cl)
        print(f"  {cl}: {n} plates")

    return plate_to_cellline


def download_plate_profiles(plate_barcode):
    """Download normalized+feature-selected profiles for one plate."""
    ts = lambda: datetime.now().strftime("%H:%M:%S")
    cache_file = CPG0000_CACHE / f"{plate_barcode}_profiles.parquet"

    if cache_file.exists():
        return pd.read_parquet(cache_file)

    url = f"{S3_BASE}/{plate_barcode}/{plate_barcode}{PROFILE_SUFFIX}"
    local_gz = CPG0000_CACHE / f"{plate_barcode}.csv.gz"

    print(f"[{ts()}] Downloading {plate_barcode}...")
    urllib.request.urlretrieve(url, local_gz)

    df = pd.read_csv(local_gz, compression="gzip")
    df.to_parquet(cache_file)
    local_gz.unlink()

    return df


def load_depmap_essentiality():
    chronos = pd.read_csv(DEPMAP_PATH, index_col=0)
    gene_breadth = {}
    n_cl = chronos.shape[0]
    for col in chronos.columns:
        gene = col.split(" ")[0]
        gene_breadth[gene] = (chronos[col] < CHRONOS_THRESHOLD).sum() / n_cl
    return pd.Series(gene_breadth, name="essentiality_breadth")


def main():
    ts = lambda: datetime.now().strftime("%H:%M:%S")
    print(f"[{ts()}] Experiment G: cpg0000 pilot cross-cell-line CRISPR")

    plate_to_cellline = get_crispr_plate_cell_lines()

    all_profiles = {"A549": [], "U2OS": []}
    for plate, cl in sorted(plate_to_cellline.items()):
        df = download_plate_profiles(plate)
        df["cell_line"] = cl
        all_profiles[cl].append(df)

    a549 = pd.concat(all_profiles["A549"], ignore_index=True)
    u2os = pd.concat(all_profiles["U2OS"], ignore_index=True)
    print(f"[{ts()}] A549: {a549.shape}, U2OS: {u2os.shape}")

    feature_cols = [
        c for c in a549.columns
        if c.startswith("Cells_") or c.startswith("Cytoplasm_") or c.startswith("Nuclei_")
    ]
    shared_features = sorted(set(feature_cols) & set(u2os.columns))
    print(f"[{ts()}] Shared features: {len(shared_features)}")

    gene_col = "Metadata_gene"
    a549_trt = a549[a549["Metadata_pert_type"] == "trt"]
    u2os_trt = u2os[u2os["Metadata_pert_type"] == "trt"]

    a549_mean = a549_trt.groupby(gene_col)[shared_features].mean()
    u2os_mean = u2os_trt.groupby(gene_col)[shared_features].mean()

    shared_genes = sorted(set(a549_mean.index) & set(u2os_mean.index))
    print(f"[{ts()}] Shared genes (treatments only): {len(shared_genes)}")

    di_values = {}
    for gene in shared_genes:
        s_a549 = a549_mean.loc[gene].values.astype(np.float64)
        s_u2os = u2os_mean.loc[gene].values.astype(np.float64)

        norm_a = np.linalg.norm(s_a549)
        norm_u = np.linalg.norm(s_u2os)
        if norm_a < 1e-10 or norm_u < 1e-10:
            continue

        cos_sim = np.dot(s_a549, s_u2os) / (norm_a * norm_u)
        cos_sim = np.clip(cos_sim, -1.0, 1.0)
        di_values[gene] = 1.0 - cos_sim

    print(f"[{ts()}] DI computed for {len(di_values)} genes (K=2, degenerate)")

    ess = load_depmap_essentiality()

    di_series = pd.Series(di_values, name="di_cross_cellline")
    merged = pd.DataFrame({"di": di_series, "essentiality_breadth": ess}).dropna()
    print(f"[{ts()}] Overlap with DepMap: {len(merged)} genes")

    if len(merged) < 20:
        result = {
            "experiment": "G_cpg0000_pilot",
            "status": "INSUFFICIENT_OVERLAP",
            "n_di_genes": len(di_values),
            "n_depmap_genes": len(ess),
            "n_overlap": len(merged),
        }
        save_results(result, RESULTS_PATH)
        return

    rho, p = stats.spearmanr(merged["di"], merged["essentiality_breadth"])
    n = len(merged)

    se = 1.0 / np.sqrt(n - 3)
    z = np.arctanh(rho)
    ci_low = float(np.tanh(z - 1.96 * se))
    ci_high = float(np.tanh(z + 1.96 * se))

    if rho < 0 and p < ALPHA:
        verdict = "directional replication"
    elif rho < 0:
        verdict = "consistent direction, underpowered"
    else:
        verdict = "does not replicate"

    print(f"[{ts()}] Spearman rho = {rho:.4f} (p = {p:.4e})")
    print(f"[{ts()}] CI = [{ci_low:.4f}, {ci_high:.4f}]")
    print(f"[{ts()}] Verdict: {verdict}")

    n_a549_wells = len(a549_trt)
    n_u2os_wells = len(u2os_trt)

    result = {
        "experiment": "G_cpg0000_pilot",
        "status": "COMPLETED",
        "n": n,
        "n_cell_lines": 2,
        "cell_lines": CELL_LINES,
        "n_a549_plates": sum(1 for v in plate_to_cellline.values() if v == "A549"),
        "n_u2os_plates": sum(1 for v in plate_to_cellline.values() if v == "U2OS"),
        "n_a549_wells": n_a549_wells,
        "n_u2os_wells": n_u2os_wells,
        "n_features": len(shared_features),
        "di_definition": "1 - cos(s_A549, s_U2OS) — degenerate K=2 case",
        "spearman_rho": round(rho, 4),
        "p_value": p,
        "ci_low": round(ci_low, 4),
        "ci_high": round(ci_high, 4),
        "alpha": ALPHA,
        "verdict": verdict,
        "batch1_comparison": {
            "batch1_rho": -0.332,
            "batch1_n": 7653,
            "batch1_context": "plate replicates within U2OS",
            "batch3_context": "across A549 and U2OS cell lines",
        },
        "limitations": [
            "K=2 cell lines: DI is a single cosine distance, not averaged over pairs",
            "No magnitude correction (n too small for meaningful OLS)",
            "No bootstrap CIs (K=2 has only one bootstrap pattern)",
            "A549 plates vary in timepoint (96h vs 144h) and antibiotics (some with Puromycin)",
            "Feature space may differ from cpg0016 production pipeline",
        ],
    }

    save_results(result, RESULTS_PATH)


if __name__ == "__main__":
    main()
