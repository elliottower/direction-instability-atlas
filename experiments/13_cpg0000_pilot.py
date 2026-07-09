"""Experiment G: cpg0000 pilot cross-cell-line CRISPR essentiality.

Tests whether the essentiality sign-flip (Batch 1 Exp C: rho=-0.33)
replicates in direction when DI is computed across cell lines (A549
and U2OS) rather than plate replicates.

cpg0000 is the JUMP pilot dataset with CRISPR profiles in two cell
lines. With K=2, DI degenerates to 1 - cos(s_A549, s_U2OS). This is
a directional convergent check, not a confirmatory test.

Pre-registered in PREREGISTRATION_BATCH3.md (Batch 3, reviewer-prompted).

Inputs:
  - cpg0000 pilot: Cell Painting Gallery (AWS S3)
  - DepMap Chronos 24Q4: data/depmap/24q4/CRISPRGeneEffect.csv

Output: results/cpg0000_pilot.json
"""
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

S3_BASE = "https://cellpainting-gallery.s3.amazonaws.com/cpg0000-jump-pilot/source_all/workspace/profiles"


def download_pilot_profiles():
    """Download cpg0000 pilot CRISPR well-level profiles."""
    CPG0000_CACHE.mkdir(parents=True, exist_ok=True)

    ts = lambda: datetime.now().strftime("%H:%M:%S")

    profiles_by_cell_line = {}

    index_url = f"{S3_BASE}/profile_index.csv"
    try:
        print(f"[{ts()}] Trying profile index at {index_url}")
        urllib.request.urlretrieve(index_url, CPG0000_CACHE / "profile_index.csv")
        index = pd.read_csv(CPG0000_CACHE / "profile_index.csv")
        print(f"[{ts()}] Profile index: {len(index)} entries")
        print(f"[{ts()}] Columns: {list(index.columns)}")
        return index
    except Exception as e:
        print(f"[{ts()}] No profile index found ({e}), trying direct plate downloads")

    for cell_line in CELL_LINES:
        cache_file = CPG0000_CACHE / f"{cell_line}_crispr_profiles.parquet"
        if cache_file.exists():
            print(f"[{ts()}] Loading cached {cell_line} profiles")
            profiles_by_cell_line[cell_line] = pd.read_parquet(cache_file)
            continue

        plate_urls = [
            f"{S3_BASE}/{cell_line}/CRISPR/profiles.parquet",
            f"{S3_BASE}/{cell_line}_CRISPR_profiles.parquet",
        ]

        for url in plate_urls:
            try:
                print(f"[{ts()}] Trying {url}")
                local_path = CPG0000_CACHE / f"{cell_line}_try.parquet"
                urllib.request.urlretrieve(url, local_path)
                df = pd.read_parquet(local_path)
                df.to_parquet(cache_file)
                profiles_by_cell_line[cell_line] = df
                print(f"[{ts()}] {cell_line}: {df.shape}")
                break
            except Exception as e:
                print(f"[{ts()}] Failed: {e}")
                continue

    return profiles_by_cell_line


def try_s3_listing():
    """List available files in cpg0000 profiles directory."""
    ts = lambda: datetime.now().strftime("%H:%M:%S")

    list_url = "https://cellpainting-gallery.s3.amazonaws.com/?prefix=cpg0000-jump-pilot/source_all/workspace/profiles/&delimiter=/&max-keys=100"
    try:
        print(f"[{ts()}] Listing S3 prefix...")
        req = urllib.request.Request(list_url)
        with urllib.request.urlopen(req) as resp:
            content = resp.read().decode("utf-8")

        prefixes = []
        for line in content.split("<CommonPrefixes>"):
            if "<Prefix>" in line:
                prefix = line.split("<Prefix>")[1].split("</Prefix>")[0]
                prefixes.append(prefix)

        keys = []
        for line in content.split("<Key>"):
            if "</Key>" in line:
                key = line.split("</Key>")[0]
                keys.append(key)

        print(f"[{ts()}] Found {len(prefixes)} subdirectories, {len(keys)} files")
        for p in prefixes[:20]:
            print(f"  {p}")
        for k in keys[:20]:
            print(f"  {k}")
        return prefixes, keys
    except Exception as e:
        print(f"[{ts()}] S3 listing failed: {e}")
        return [], []


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

    prefixes, keys = try_s3_listing()

    if not prefixes and not keys:
        print(f"[{ts()}] Cannot access cpg0000 data. Checking if cached profiles exist...")

    cached_a549 = CPG0000_CACHE / "A549_crispr_profiles.parquet"
    cached_u2os = CPG0000_CACHE / "U2OS_crispr_profiles.parquet"

    if cached_a549.exists() and cached_u2os.exists():
        a549 = pd.read_parquet(cached_a549)
        u2os = pd.read_parquet(cached_u2os)
        print(f"[{ts()}] Loaded cached: A549 {a549.shape}, U2OS {u2os.shape}")
    else:
        profiles = download_pilot_profiles()
        if isinstance(profiles, dict) and len(profiles) == 2:
            a549 = profiles["A549"]
            u2os = profiles["U2OS"]
        else:
            print(f"[{ts()}] Could not load cpg0000 profiles. S3 structure:")
            for p in prefixes:
                print(f"  {p}")
            for k in keys[:30]:
                print(f"  {k}")

            result = {
                "experiment": "G_cpg0000_pilot",
                "status": "DATA_UNAVAILABLE",
                "note": "cpg0000 pilot profiles could not be downloaded. S3 structure may require aws cli or different path format.",
                "s3_prefixes_found": prefixes[:20],
                "s3_keys_found": keys[:30],
            }
            save_results(result, RESULTS_PATH)
            return

    feature_cols_a549 = [c for c in a549.columns if c.startswith("X_") or c.startswith("PC_") or c.startswith("Cells_") or c.startswith("Cytoplasm_") or c.startswith("Nuclei_")]
    feature_cols_u2os = [c for c in u2os.columns if c.startswith("X_") or c.startswith("PC_") or c.startswith("Cells_") or c.startswith("Cytoplasm_") or c.startswith("Nuclei_")]

    shared_features = sorted(set(feature_cols_a549) & set(feature_cols_u2os))
    print(f"[{ts()}] Shared features: {len(shared_features)}")

    if not shared_features:
        print(f"[{ts()}] No shared feature columns found")
        print(f"  A549 columns sample: {list(a549.columns[:10])}")
        print(f"  U2OS columns sample: {list(u2os.columns[:10])}")
        result = {
            "experiment": "G_cpg0000_pilot",
            "status": "NO_SHARED_FEATURES",
            "a549_columns_sample": list(a549.columns[:20]),
            "u2os_columns_sample": list(u2os.columns[:20]),
        }
        save_results(result, RESULTS_PATH)
        return

    gene_col = None
    for candidate in ["Metadata_Symbol", "Metadata_gene", "Metadata_target", "Metadata_JCP2022"]:
        if candidate in a549.columns and candidate in u2os.columns:
            gene_col = candidate
            break

    if gene_col is None:
        print(f"[{ts()}] No gene column found in pilot data")
        result = {
            "experiment": "G_cpg0000_pilot",
            "status": "NO_GENE_COLUMN",
            "a549_meta_columns": [c for c in a549.columns if "Meta" in c or "gene" in c.lower()],
            "u2os_meta_columns": [c for c in u2os.columns if "Meta" in c or "gene" in c.lower()],
        }
        save_results(result, RESULTS_PATH)
        return

    a549_mean = a549.groupby(gene_col)[shared_features].mean()
    u2os_mean = u2os.groupby(gene_col)[shared_features].mean()

    shared_genes = sorted(set(a549_mean.index) & set(u2os_mean.index))
    print(f"[{ts()}] Shared genes: {len(shared_genes)}")

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

    result = {
        "experiment": "G_cpg0000_pilot",
        "status": "COMPLETED",
        "n": n,
        "n_cell_lines": 2,
        "cell_lines": CELL_LINES,
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
            "Feature space may differ from cpg0016 production",
        ],
    }

    save_results(result, RESULTS_PATH)


if __name__ == "__main__":
    main()
