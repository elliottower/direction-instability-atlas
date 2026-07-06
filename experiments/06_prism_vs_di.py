"""Experiment PRISM (EXPLORATORY): PRISM viability CV vs LINCS DI.

Tests whether drugs with low DI in LINCS show low viability CV in PRISM.
Spearman rho > 0.15 for shared drugs (n~1,726).

EXPLORATORY — uncorrected p-value, not a primary test.

Inputs:
  - LINCS DI: drug-perturbation-geometry/zenodo_v1/drug_instability_8949.csv
  - PRISM CV: results/prism_cv.csv (from script 04)

Output: results/prism_concordance.json
"""
from datetime import datetime
from pathlib import Path

import pandas as pd
from scipy import stats

from experiments.utils import ci_lower_bound, ci_upper_bound, save_results

LINCS_PATH = Path("../drug-perturbation-geometry/zenodo_v1/drug_instability_8949.csv")
PRISM_CV_PATH = Path("results/prism_cv.csv")
RESULTS_PATH = Path("results/prism_concordance.json")


def main():
    ts = lambda: datetime.now().strftime("%H:%M:%S")
    print(f"[{ts()}] PRISM exploratory: viability CV vs LINCS DI")

    lincs = pd.read_csv(LINCS_PATH)
    prism = pd.read_csv(PRISM_CV_PATH)
    print(f"[{ts()}] LINCS: {len(lincs)} drugs, PRISM: {len(prism)} drugs")

    lincs["drug_lower"] = lincs["drug_name"].str.lower().str.strip()
    prism["drug_lower"] = prism["drug_name"].str.lower().str.strip()

    merged = lincs.merge(prism, on="drug_lower", suffixes=("_lincs", "_prism"))
    print(f"[{ts()}] Overlap: {len(merged)} drugs")

    if len(merged) == 0:
        save_results({"experiment": "PRISM_exploratory", "n": 0,
                       "note": "No drug overlap found"}, RESULTS_PATH)
        return

    rho, p = stats.spearmanr(merged["direction_instability"], merged["viability_cv"])
    n = len(merged)

    result = {
        "experiment": "PRISM_exploratory",
        "label": "EXPLORATORY — uncorrected p-value",
        "n": n,
        "spearman_rho": round(rho, 4),
        "p_value": p,
        "ci_low": round(ci_lower_bound(rho, n), 4),
        "ci_high": round(ci_upper_bound(rho, n), 4),
        "threshold": 0.15,
        "passes_threshold": bool(rho > 0.15),
    }

    print(f"[{ts()}] rho = {rho:.4f}, p = {p:.2e}, n = {n}")
    print(f"[{ts()}] {'PASS' if rho > 0.15 else 'FAIL'} (threshold: rho > 0.15)")

    save_results(result, RESULTS_PATH)


if __name__ == "__main__":
    main()
