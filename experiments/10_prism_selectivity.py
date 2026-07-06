"""Experiment E3: Drug selectivity — does DI predict selective killing?

Tests whether high-DI drugs are selective killers (Gini of kill
pattern) vs. low-DI drugs being broad/uniform killers.

Primary metric: gini_corrected (Gini residualized against mean kill
magnitude via OLS — same pattern as di_corrected). Raw Gini on kill
depth is a magnitude-variance confound trap.

Partial Spearman rho (Pearson-residualize then Spearman) controlling
for LINCS n_contexts (1 covariate). Mean kill magnitude is already
removed from Gini via the correction, NOT double-counted as a
covariate.
Decision criterion: CI lower bound > 0.05, p < 0.0167.

Inputs:
  - LINCS DI: drug-perturbation-geometry/zenodo_v1/drug_instability_8949.csv
  - PRISM LFC: data/prism/Repurposing_Public_24Q2_LFC_COLLAPSED.csv
  - PRISM compounds: data/prism/Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv

Output: results/prism_selectivity.json
"""
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LinearRegression

from experiments.utils import (
    ci_lower_bound,
    ci_upper_bound,
    fisher_z_pvalue,
    gini_coefficient,
    partial_spearman,
    save_results,
)

LINCS_PATH = Path("../drug-perturbation-geometry/zenodo_v1/drug_instability_8949.csv")
LFC_PATH = Path("data/prism/Repurposing_Public_24Q2_LFC_COLLAPSED.csv")
COMPOUND_PATH = Path("data/prism/Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv")
RESULTS_PATH = Path("results/prism_selectivity.json")

ALPHA = 0.0167
CI_FLOOR = 0.05
MIN_CELL_LINES = 20
KILL_THRESHOLD = -0.5


def compute_drug_selectivity(lfc_df: pd.DataFrame, compound_meta: pd.DataFrame) -> pd.DataFrame:
    compound_meta = compound_meta.copy()
    compound_meta["broad_id"] = compound_meta["IDs"].str.replace("BRD:", "", regex=False)
    brd_to_name = dict(zip(compound_meta["broad_id"], compound_meta["Drug.Name"]))

    lfc_df = lfc_df[~lfc_df["broad_id"].str.contains("QC Failure", na=False)]

    rows = []
    for broad_id, group in lfc_df.groupby("broad_id"):
        mean_lfc_per_cell = group.groupby("row_id")["LFC"].mean().dropna()
        lfc_values = mean_lfc_per_cell.values

        if len(lfc_values) < MIN_CELL_LINES:
            continue

        mean_lfc = float(np.mean(lfc_values))
        if abs(mean_lfc) < 0.01:
            continue

        neg_lfc = lfc_values[lfc_values < 0]
        if len(neg_lfc) < 5:
            continue

        gini_raw = gini_coefficient(neg_lfc)
        mean_kill_magnitude = float(np.mean(np.abs(neg_lfc)))
        kill_fraction = float(np.mean(lfc_values < KILL_THRESHOLD))

        drug_name = brd_to_name.get(broad_id, broad_id)

        rows.append({
            "broad_id": broad_id,
            "drug_name": drug_name,
            "gini_raw": gini_raw,
            "mean_kill_magnitude": mean_kill_magnitude,
            "kill_fraction": kill_fraction,
            "mean_lfc": mean_lfc,
            "n_cell_lines": len(lfc_values),
        })

    df = pd.DataFrame(rows)

    if len(df) > 0:
        X = df["mean_kill_magnitude"].values.reshape(-1, 1)
        y = df["gini_raw"].values
        reg = LinearRegression().fit(X, y)
        df["gini_corrected"] = y - reg.predict(X) + reg.intercept_

    return df


def main():
    ts = lambda: datetime.now().strftime("%H:%M:%S")
    print(f"[{ts()}] Experiment E3: Drug selectivity (Gini) vs LINCS DI")

    lincs = pd.read_csv(LINCS_PATH)
    print(f"[{ts()}] LINCS DI: {len(lincs)} drugs")

    lfc_df = pd.read_csv(LFC_PATH)
    compound_meta = pd.read_csv(COMPOUND_PATH)
    print(f"[{ts()}] PRISM LFC loaded: {lfc_df.shape[0]} rows")

    selectivity = compute_drug_selectivity(lfc_df, compound_meta)
    print(f"[{ts()}] Selectivity computed: {len(selectivity)} drugs")

    lincs["drug_lower"] = lincs["drug_name"].str.lower().str.strip()
    selectivity["drug_lower"] = selectivity["drug_name"].str.lower().str.strip()

    merged = lincs.merge(selectivity, on="drug_lower", suffixes=("_lincs", "_prism"))
    print(f"[{ts()}] Overlap: {len(merged)} drugs")

    if len(merged) == 0:
        save_results({"experiment": "E3_prism_selectivity", "n": 0}, RESULTS_PATH)
        return

    di_vals = merged["direction_instability"].values.astype(np.float64)
    gini_corr = merged["gini_corrected"].values.astype(np.float64)
    gini_raw_vals = merged["gini_raw"].values.astype(np.float64)
    kill_frac = merged["kill_fraction"].values.astype(np.float64)
    n_ctx = merged["n_cell_lines_lincs"].values.astype(np.float64).reshape(-1, 1)
    n = len(merged)

    # Primary: gini_corrected
    raw_rho, raw_p = stats.spearmanr(di_vals, gini_corr)
    partial_rho = partial_spearman(di_vals, gini_corr, n_ctx)
    ci_low = ci_lower_bound(partial_rho, n, n_cov=1)
    ci_high = ci_upper_bound(partial_rho, n, n_cov=1)
    p_val = fisher_z_pvalue(partial_rho, n, n_cov=1)
    passes = ci_low > CI_FLOOR and p_val < ALPHA

    # Sensitivity: gini_raw (NOT primary — magnitude-confounded)
    raw_rho_uncorr, _ = stats.spearmanr(di_vals, gini_raw_vals)
    partial_rho_uncorr = partial_spearman(di_vals, gini_raw_vals, n_ctx)

    # Sensitivity: kill_fraction
    raw_rho_kf, _ = stats.spearmanr(di_vals, kill_frac)
    partial_rho_kf = partial_spearman(di_vals, kill_frac, n_ctx)

    if abs(partial_rho) >= 0.20:
        interpretation = "strong: DI predicts drug selectivity"
    elif abs(partial_rho) >= 0.05:
        interpretation = "moderate: some selectivity signal"
    else:
        interpretation = "null: DI does not predict selective vs. broad killing"

    result = {
        "experiment": "E3_prism_selectivity",
        "n": n,
        "primary_metric": "gini_corrected",
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
        "magnitude_correction": "Gini residualized against mean_kill_magnitude via OLS",
        "sensitivity_gini_raw": {
            "raw_rho": round(raw_rho_uncorr, 4),
            "partial_rho": round(partial_rho_uncorr, 4),
            "note": "NOT primary — magnitude-confounded",
        },
        "sensitivity_kill_fraction": {
            "raw_rho": round(raw_rho_kf, 4),
            "partial_rho": round(partial_rho_kf, 4),
            "threshold": KILL_THRESHOLD,
        },
    }

    print(f"[{ts()}] Primary (gini_corrected): partial rho = {partial_rho:.4f}")
    print(f"[{ts()}] CI = [{ci_low:.4f}, {ci_high:.4f}], p = {p_val:.2e}")
    print(f"[{ts()}] Sensitivity (gini_raw): partial rho = {partial_rho_uncorr:.4f}")
    print(f"[{ts()}] Sensitivity (kill_fraction): partial rho = {partial_rho_kf:.4f}")
    print(f"[{ts()}] {'PASS' if passes else 'FAIL'} — {interpretation}")

    save_results(result, RESULTS_PATH)


if __name__ == "__main__":
    main()
