"""Experiment D: Knockout vs overexpression DI correlation.

Tests whether genes with high CRISPR knockout DI also have high ORF
overexpression DI. Context-dependence should be a property of the gene,
not the perturbation direction.

Universe: 5,220 genes with >= 5 replicates in both CRISPR and ORF.

Decision criterion: CI lower > 0.05 AND observed rho exceeds 95th
percentile of permutation null (1,000 shuffles). Alpha = 0.0125.

Includes exploratory directional divergence analysis.

Inputs:
  - CRISPR DI: results/crispr_di.csv (from script 02)
  - ORF DI: results/orf_di.csv (from script 03)

Output: results/knockout_vs_overexpression.json
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
    save_results,
)

CRISPR_DI_PATH = Path("results/crispr_di.csv")
ORF_DI_PATH = Path("results/orf_di.csv")
RESULTS_PATH = Path("results/knockout_vs_overexpression.json")

ALPHA = 0.0125
CI_FLOOR = 0.05
N_PERMUTATIONS = 1000


def permutation_null(
    crispr_di: np.ndarray,
    orf_di: np.ndarray,
    n_perms: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Shuffle ORF gene labels and recompute Spearman rho n_perms times."""
    null_rhos = np.empty(n_perms)
    for i in range(n_perms):
        shuffled = rng.permutation(orf_di)
        null_rhos[i], _ = stats.spearmanr(crispr_di, shuffled)
    return null_rhos


def directional_divergence(merged: pd.DataFrame) -> dict:
    """Exploratory: classify genes where |DI_CRISPR - DI_ORF| > 1 SD."""
    di_crispr = merged["di_corrected_crispr"].values
    di_orf = merged["di_corrected_orf"].values
    diff = di_crispr - di_orf
    sd = np.std(diff, ddof=1)

    ko_stable_oe_unstable = merged[diff < -sd]["gene"].tolist()
    ko_unstable_oe_stable = merged[diff > sd]["gene"].tolist()

    return {
        "label": "EXPLORATORY — no pass/fail threshold",
        "sd_of_difference": round(float(sd), 4),
        "knockout_stable_overexpression_unstable": {
            "count": len(ko_stable_oe_unstable),
            "predicted_biology": "dosage-sensitive gain-of-function genes",
            "genes": ko_stable_oe_unstable[:20],
        },
        "knockout_unstable_overexpression_stable": {
            "count": len(ko_unstable_oe_stable),
            "predicted_biology": "genes with redundant paralogs",
            "genes": ko_unstable_oe_stable[:20],
        },
    }


def main():
    ts = lambda: datetime.now().strftime("%H:%M:%S")
    print(f"[{ts()}] Experiment D: Knockout vs overexpression DI")

    crispr = pd.read_csv(CRISPR_DI_PATH)
    orf = pd.read_csv(ORF_DI_PATH)
    print(f"[{ts()}] CRISPR: {len(crispr)} genes, ORF: {len(orf)} genes")

    merged = crispr.merge(orf, on="gene", suffixes=("_crispr", "_orf"))
    print(f"[{ts()}] Overlap: {len(merged)} genes")

    if len(merged) == 0:
        save_results({"experiment": "D_knockout_vs_overexpression", "n": 0}, RESULTS_PATH)
        return

    di_crispr = merged["di_corrected_crispr"].values.astype(np.float64)
    di_orf = merged["di_corrected_orf"].values.astype(np.float64)
    n = len(merged)

    observed_rho, observed_p = stats.spearmanr(di_crispr, di_orf)
    ci_low = ci_lower_bound(observed_rho, n)
    ci_high = ci_upper_bound(observed_rho, n)
    p_val = fisher_z_pvalue(observed_rho, n)

    # Permutation null
    print(f"[{ts()}] Running {N_PERMUTATIONS} permutations")
    rng = np.random.default_rng(42)
    null_rhos = permutation_null(di_crispr, di_orf, N_PERMUTATIONS, rng)

    perm_p = float(np.mean(null_rhos >= observed_rho))
    perm_95th = float(np.percentile(null_rhos, 95))

    passes_ci = bool(ci_low > CI_FLOOR)
    passes_perm = bool(observed_rho > perm_95th)
    passes = passes_ci and passes_perm and p_val < ALPHA

    print(f"[{ts()}] Observed rho = {observed_rho:.4f}")
    print(f"[{ts()}] CI = [{ci_low:.4f}, {ci_high:.4f}]")
    print(f"[{ts()}] Permutation null: mean={null_rhos.mean():.4f}, "
          f"SD={null_rhos.std():.4f}, 95th={perm_95th:.4f}")
    print(f"[{ts()}] CI lower > {CI_FLOOR}: {'PASS' if passes_ci else 'FAIL'}")
    print(f"[{ts()}] Exceeds permutation 95th: {'PASS' if passes_perm else 'FAIL'}")
    print(f"[{ts()}] Overall: {'PASS' if passes else 'FAIL'}")

    # Exploratory: directional divergence
    print(f"[{ts()}] Directional divergence analysis (exploratory)")
    divergence = directional_divergence(merged)

    result = {
        "experiment": "D_knockout_vs_overexpression",
        "n": n,
        "observed_spearman_rho": round(observed_rho, 4),
        "ci_low": round(ci_low, 4),
        "ci_high": round(ci_high, 4),
        "fisher_z_p_value": p_val,
        "permutation_null": {
            "n_permutations": N_PERMUTATIONS,
            "mean": round(float(null_rhos.mean()), 4),
            "sd": round(float(null_rhos.std()), 4),
            "percentile_95th": round(perm_95th, 4),
            "empirical_p_value": perm_p,
        },
        "criterion_1": f"CI lower > {CI_FLOOR}",
        "criterion_2": "observed rho > permutation 95th percentile",
        "alpha": ALPHA,
        "passes_ci": passes_ci,
        "passes_permutation": passes_perm,
        "passes_overall": passes,
        "directional_divergence": divergence,
    }

    save_results(result, RESULTS_PATH)


if __name__ == "__main__":
    main()
