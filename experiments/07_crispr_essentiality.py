"""Experiment C: CRISPR DI predicts essentiality.

Tests whether genes essential in more cell types (DepMap Chronos < -0.5)
show higher morphological DI in JUMP-CP CRISPR knockouts.

Partial Spearman rho controlling for n_contexts (1 covariate).
Decision criterion: CI lower bound > 0.05, p < 0.0125.

Inputs:
  - CRISPR DI: results/crispr_di.csv (from script 02)
  - DepMap Chronos: downloaded from DepMap portal

Output: results/essentiality.json
"""
import urllib.request
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from experiments.utils import (
    ci_lower_bound,
    ci_upper_bound,
    fisher_z_pvalue,
    partial_spearman,
    save_results,
)

CRISPR_DI_PATH = Path("results/crispr_di.csv")
DEPMAP_URL = "https://ndownloader.figshare.com/files/51064667"
DEPMAP_CACHE = Path("data/depmap/24q4/CRISPRGeneEffect.csv")
RESULTS_PATH = Path("results/essentiality.json")

ALPHA = 0.0125
CI_FLOOR = 0.05
CHRONOS_THRESHOLD = -0.5


def load_depmap_essentiality() -> pd.DataFrame:
    """Load DepMap Chronos scores and compute essentiality breadth per gene."""
    if DEPMAP_CACHE.exists():
        chronos = pd.read_csv(DEPMAP_CACHE, index_col=0)
    else:
        DEPMAP_CACHE.parent.mkdir(parents=True, exist_ok=True)
        print(f"Downloading DepMap Chronos to {DEPMAP_CACHE}")
        urllib.request.urlretrieve(DEPMAP_URL, DEPMAP_CACHE)
        chronos = pd.read_csv(DEPMAP_CACHE, index_col=0)

    # Columns are "GENE (ENTREZ_ID)" format — extract gene symbol
    gene_breadth = {}
    n_cell_lines = chronos.shape[0]

    for col in chronos.columns:
        gene = col.split(" ")[0]
        essential_count = (chronos[col] < CHRONOS_THRESHOLD).sum()
        gene_breadth[gene] = essential_count / n_cell_lines

    return pd.DataFrame([
        {"gene": gene, "essentiality_breadth": breadth}
        for gene, breadth in gene_breadth.items()
    ])


def main():
    ts = lambda: datetime.now().strftime("%H:%M:%S")
    print(f"[{ts()}] Experiment C: CRISPR DI vs essentiality")

    crispr_di = pd.read_csv(CRISPR_DI_PATH)
    print(f"[{ts()}] CRISPR DI: {len(crispr_di)} genes")

    essentiality = load_depmap_essentiality()
    print(f"[{ts()}] DepMap essentiality: {len(essentiality)} genes")

    merged = crispr_di.merge(essentiality, on="gene")
    print(f"[{ts()}] Overlap: {len(merged)} genes")

    if len(merged) == 0:
        save_results({"experiment": "C_essentiality", "n": 0}, RESULTS_PATH)
        return

    di_vals = merged["di_corrected"].values.astype(np.float64)
    ess_vals = merged["essentiality_breadth"].values.astype(np.float64)
    n_ctx = merged["n_contexts"].values.astype(np.float64).reshape(-1, 1)
    n = len(merged)

    raw_rho, raw_p = stats.spearmanr(di_vals, ess_vals)
    partial_rho = partial_spearman(di_vals, ess_vals, n_ctx)

    ci_low = ci_lower_bound(partial_rho, n, n_cov=1)
    ci_high = ci_upper_bound(partial_rho, n, n_cov=1)
    p_val = fisher_z_pvalue(partial_rho, n, n_cov=1)

    passes = ci_low > CI_FLOOR and p_val < ALPHA

    if partial_rho >= 0.15:
        interpretation = "biologically meaningful"
    elif partial_rho >= 0.05:
        interpretation = "weak — statistically detectable but explaining < 2% of variance"
    else:
        interpretation = "negligible"

    result = {
        "experiment": "C_essentiality",
        "n": n,
        "raw_spearman_rho": round(raw_rho, 4),
        "raw_p_value": raw_p,
        "partial_spearman_rho": round(partial_rho, 4),
        "partial_ci_low": round(ci_low, 4),
        "partial_ci_high": round(ci_high, 4),
        "partial_p_value": p_val,
        "criterion": f"CI lower > {CI_FLOOR}",
        "alpha": ALPHA,
        "passes": passes,
        "interpretation": interpretation,
        "method": "Pearson-residualize then Spearman, 1 covariate (n_contexts)",
        "essentiality_definition": f"fraction of cell lines with Chronos < {CHRONOS_THRESHOLD}",
    }

    print(f"[{ts()}] Raw rho = {raw_rho:.4f}, Partial rho = {partial_rho:.4f}")
    print(f"[{ts()}] CI = [{ci_low:.4f}, {ci_high:.4f}], p = {p_val:.2e}")
    print(f"[{ts()}] {'PASS' if passes else 'FAIL'} — {interpretation}")

    save_results(result, RESULTS_PATH)


if __name__ == "__main__":
    main()
