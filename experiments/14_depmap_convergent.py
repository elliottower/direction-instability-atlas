"""Experiment H: DepMap essentiality variability as convergent evidence.

Tests whether broadly essential genes show lower cross-cell-line
variability in Chronos dependency scores than selectively essential
genes. Convergent evidence for the biological hypothesis (functional
constraint suppresses context-dependence), using a DIFFERENT geometric
construct from morphological DI.

Pre-registered in PREREGISTRATION_BATCH3.md (Batch 3, reviewer-prompted).

Inputs:
  - DepMap Chronos 24Q4: data/depmap/24q4/CRISPRGeneEffect.csv

Output: results/depmap_convergent.json
"""
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

DEPMAP_PATH = Path("data/depmap/24q4/CRISPRGeneEffect.csv")
RESULTS_PATH = Path("results/depmap_convergent.json")

ALPHA = 0.0167
CI_FLOOR = 0.05
CHRONOS_THRESHOLD = -0.5


def main():
    ts = lambda: datetime.now().strftime("%H:%M:%S")
    print(f"[{ts()}] Experiment H: DepMap essentiality variability (convergent)")

    chronos = pd.read_csv(DEPMAP_PATH, index_col=0)
    n_cell_lines = chronos.shape[0]
    print(f"[{ts()}] DepMap Chronos: {chronos.shape[0]} cell lines, {chronos.shape[1]} genes")

    rows = []
    for col in chronos.columns:
        gene = col.split(" ")[0]
        scores = chronos[col].dropna().values
        n_valid = len(scores)
        if n_valid < 50:
            continue

        breadth = (scores < CHRONOS_THRESHOLD).sum() / n_valid
        variability = float(np.std(scores))

        rows.append({
            "gene": gene,
            "essentiality_breadth": breadth,
            "dependency_sd": variability,
            "n_valid_cell_lines": n_valid,
        })

    df = pd.DataFrame(rows)
    print(f"[{ts()}] {len(df)} genes with >= 50 valid Chronos values")

    breadth = df["essentiality_breadth"].values.astype(np.float64)
    variability = df["dependency_sd"].values.astype(np.float64)
    n_valid = df["n_valid_cell_lines"].values.astype(np.float64).reshape(-1, 1)

    raw_rho, raw_p = stats.spearmanr(variability, breadth)

    partial_rho = partial_spearman(variability, breadth, n_valid)
    n = len(df)

    ci_low = ci_lower_bound(partial_rho, n, n_cov=1)
    ci_high = ci_upper_bound(partial_rho, n, n_cov=1)
    p_val = fisher_z_pvalue(partial_rho, n, n_cov=1)

    passes = ci_high < -CI_FLOOR and p_val < ALPHA

    print(f"[{ts()}] Raw Spearman rho = {raw_rho:.4f} (p = {raw_p:.2e})")
    print(f"[{ts()}] Partial Spearman rho = {partial_rho:.4f}")
    print(f"[{ts()}] CI = [{ci_low:.4f}, {ci_high:.4f}], p = {p_val:.2e}")
    print(f"[{ts()}] {'PASS' if passes else 'FAIL'}")

    if partial_rho < -0.15:
        interpretation = "convergent: functional constraint suppresses context-dependence in dependency scores"
    elif partial_rho < -0.05:
        interpretation = "weakly convergent: same direction as morphological sign-flip but small effect"
    elif partial_rho > 0.05:
        interpretation = "does not converge: opposite direction from morphological sign-flip"
    else:
        interpretation = "null: no relationship between essentiality breadth and dependency variability"

    result = {
        "experiment": "H_depmap_convergent",
        "estimand": "SD of Chronos dependency scores across cell lines (DIFFERENT from morphological DI)",
        "n": n,
        "n_cell_lines_depmap": int(n_cell_lines),
        "min_valid_cell_lines": 50,
        "raw_spearman_rho": round(raw_rho, 4),
        "raw_p_value": raw_p,
        "partial_spearman_rho": round(partial_rho, 4),
        "partial_ci_low": round(ci_low, 4),
        "partial_ci_high": round(ci_high, 4),
        "partial_p_value": p_val,
        "covariate": "n_valid_cell_lines",
        "criterion": f"CI upper < -{CI_FLOOR} (predicted negative)",
        "alpha": ALPHA,
        "passes": passes,
        "interpretation": interpretation,
        "convergence_note": "This tests the biological hypothesis (functional constraint suppresses context-dependence), NOT the morphological DI metric. Different geometric construct (scalar SD vs vector cosine distance).",
    }

    save_results(result, RESULTS_PATH)


if __name__ == "__main__":
    main()
