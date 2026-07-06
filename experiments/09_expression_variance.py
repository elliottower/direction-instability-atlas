"""Experiment E2: Expression variance as Waddington basin proxy.

Tests whether genes whose expression varies more across cell lines
(shallow attractor basins) produce higher morphological DI when
knocked out.

Partial Spearman rho (Pearson-residualize then Spearman) controlling
for n_contexts (1 covariate).
Decision criterion: CI lower bound > 0.05, p < 0.0167.

Inputs:
  - CRISPR DI: results/crispr_di.csv (from script 02)
  - DepMap expression: OmicsExpressionProteinCodingGenesTPMLogp1.csv

Output: results/expression_variance.json
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
EXPRESSION_URL = "https://figshare.com/ndownloader/files/44318028"
EXPRESSION_CACHE = Path("data/depmap/OmicsExpressionProteinCodingGenesTPMLogp1.csv")
RESULTS_PATH = Path("results/expression_variance.json")

ALPHA = 0.0167
CI_FLOOR = 0.05


def load_expression_variance() -> pd.DataFrame:
    if EXPRESSION_CACHE.exists():
        expr = pd.read_csv(EXPRESSION_CACHE, index_col=0)
    else:
        EXPRESSION_CACHE.parent.mkdir(parents=True, exist_ok=True)
        print(f"Downloading DepMap expression to {EXPRESSION_CACHE}")
        urllib.request.urlretrieve(EXPRESSION_URL, EXPRESSION_CACHE)
        expr = pd.read_csv(EXPRESSION_CACHE, index_col=0)

    rows = []
    for col in expr.columns:
        gene = col.split(" ")[0]
        values = expr[col].dropna().values
        if len(values) < 10:
            continue
        rows.append({
            "gene": gene,
            "expr_variance": float(np.var(values, ddof=1)),
            "expr_mean": float(np.mean(values)),
            "expr_n_cell_lines": len(values),
        })

    return pd.DataFrame(rows)


def main():
    ts = lambda: datetime.now().strftime("%H:%M:%S")
    print(f"[{ts()}] Experiment E2: Expression variance vs CRISPR DI")

    crispr_di = pd.read_csv(CRISPR_DI_PATH)
    print(f"[{ts()}] CRISPR DI: {len(crispr_di)} genes")

    expr_var = load_expression_variance()
    print(f"[{ts()}] Expression variance: {len(expr_var)} genes")

    merged = crispr_di.merge(expr_var, on="gene")
    print(f"[{ts()}] Overlap: {len(merged)} genes")

    if len(merged) == 0:
        save_results({"experiment": "E2_expression_variance", "n": 0}, RESULTS_PATH)
        return

    di_vals = merged["di_corrected"].values.astype(np.float64)
    var_vals = merged["expr_variance"].values.astype(np.float64)
    n_ctx = merged["n_contexts"].values.astype(np.float64).reshape(-1, 1)
    n = len(merged)

    raw_rho, raw_p = stats.spearmanr(di_vals, var_vals)
    partial_rho = partial_spearman(di_vals, var_vals, n_ctx)

    ci_low = ci_lower_bound(partial_rho, n, n_cov=1)
    ci_high = ci_upper_bound(partial_rho, n, n_cov=1)
    p_val = fisher_z_pvalue(partial_rho, n, n_cov=1)

    passes = ci_low > CI_FLOOR and p_val < ALPHA

    if abs(partial_rho) >= 0.15:
        interpretation = "strong support: expression landscape predicts DI"
    elif abs(partial_rho) >= 0.05:
        interpretation = "weak but detectable"
    else:
        interpretation = "null: DI not driven by expression-level context-dependence"

    result = {
        "experiment": "E2_expression_variance",
        "n": n,
        "raw_spearman_rho": round(raw_rho, 4),
        "raw_p_value": raw_p,
        "partial_spearman_rho": round(partial_rho, 4),
        "partial_ci_low": round(ci_low, 4),
        "partial_ci_high": round(ci_high, 4),
        "partial_p_value": p_val,
        "criterion": f"CI lower > {CI_FLOOR} AND p < {ALPHA}",
        "alpha": ALPHA,
        "passes": passes,
        "interpretation": interpretation,
        "method": "Pearson-residualize then Spearman, 1 covariate (n_contexts)",
        "expression_source": "DepMap OmicsExpressionProteinCodingGenesTPMLogp1",
    }

    print(f"[{ts()}] Raw rho = {raw_rho:.4f}, Partial rho = {partial_rho:.4f}")
    print(f"[{ts()}] CI = [{ci_low:.4f}, {ci_high:.4f}], p = {p_val:.2e}")
    print(f"[{ts()}] {'PASS' if passes else 'FAIL'} — {interpretation}")

    save_results(result, RESULTS_PATH)


if __name__ == "__main__":
    main()
