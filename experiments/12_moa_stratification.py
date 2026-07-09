"""Experiment F: MoA-stratified DI characterization.

Tests whether DI varies systematically across drug mechanism-of-action
classes, after controlling for profiling depth (n_cell_lines) and
signal strength (mean_norm). Permutation null controls for class-size
imbalance.

Pre-registered in PREREGISTRATION_BATCH3.md (Batch 3, reviewer-prompted).

Inputs:
  - LINCS DI: drug-perturbation-geometry/zenodo_v1/drug_instability_8949.csv

Output: results/moa_stratification.json
"""
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LinearRegression

from experiments.utils import save_results

LINCS_DI_PATH = Path("../drug-perturbation-geometry/zenodo_v1/drug_instability_8949.csv")
RESULTS_PATH = Path("results/moa_stratification.json")
MIN_CLASS_SIZE = 5
N_PERMUTATIONS = 1000
ALPHA = 0.0167


def residualize_di(df):
    X = df[["n_cell_lines", "mean_norm"]].values
    y = df["direction_instability"].values
    reg = LinearRegression().fit(X, y)
    return y - reg.predict(X) + reg.intercept_


def compute_eta_squared(di_vals, labels):
    groups = {}
    for di, lab in zip(di_vals, labels):
        groups.setdefault(lab, []).append(di)
    group_arrays = [np.array(v) for v in groups.values()]

    grand_mean = np.mean(di_vals)
    ss_between = sum(
        len(g) * (g.mean() - grand_mean) ** 2 for g in group_arrays
    )
    ss_total = np.sum((di_vals - grand_mean) ** 2)
    if ss_total == 0:
        return 0.0
    return ss_between / ss_total


def main():
    ts = lambda: datetime.now().strftime("%H:%M:%S")
    print(f"[{ts()}] Experiment F: MoA-stratified DI characterization")

    di = pd.read_csv(LINCS_DI_PATH)
    print(f"[{ts()}] LINCS DI: {len(di)} drugs")

    moa_drugs = di[di["moa"].notna()].copy()
    print(f"[{ts()}] Drugs with MoA: {len(moa_drugs)}")

    moa_counts = moa_drugs["moa"].value_counts()
    big_moas = moa_counts[moa_counts >= MIN_CLASS_SIZE].index
    moa_sub = moa_drugs[moa_drugs["moa"].isin(big_moas)].copy()
    print(f"[{ts()}] MoAs with >= {MIN_CLASS_SIZE} drugs: {len(big_moas)} classes, {len(moa_sub)} drugs")

    moa_sub["di_residualized"] = residualize_di(moa_sub)

    di_resid = moa_sub["di_residualized"].values
    moa_labels = moa_sub["moa"].values

    raw_eta = compute_eta_squared(
        moa_sub["direction_instability"].values, moa_labels
    )
    resid_eta = compute_eta_squared(di_resid, moa_labels)

    groups_resid = [g["di_residualized"].values for _, g in moa_sub.groupby("moa")]
    f_stat, f_pval = stats.f_oneway(*groups_resid)

    print(f"[{ts()}] Raw eta^2 = {raw_eta:.4f}")
    print(f"[{ts()}] Residualized eta^2 = {resid_eta:.4f}")
    print(f"[{ts()}] F = {f_stat:.2f}, p = {f_pval:.2e}")

    print(f"[{ts()}] Permutation null ({N_PERMUTATIONS} shuffles)...")
    rng = np.random.default_rng(2026)
    null_etas = np.empty(N_PERMUTATIONS)
    for i in range(N_PERMUTATIONS):
        shuffled_labels = rng.permutation(moa_labels)
        null_etas[i] = compute_eta_squared(di_resid, shuffled_labels)

    perm_p = float((null_etas >= resid_eta).mean())
    print(f"[{ts()}] Null eta^2: mean = {null_etas.mean():.4f}, "
          f"SD = {null_etas.std():.4f}, 95th = {np.percentile(null_etas, 95):.4f}, "
          f"99th = {np.percentile(null_etas, 99):.4f}")
    print(f"[{ts()}] Permutation p = {perm_p}")

    passes = perm_p < ALPHA and resid_eta > np.percentile(null_etas, 99)

    moa_stats = moa_sub.groupby("moa").agg(
        mean_di_resid=("di_residualized", "mean"),
        std_di_resid=("di_residualized", "std"),
        mean_di_raw=("direction_instability", "mean"),
        n=("direction_instability", "count"),
    ).reset_index()
    moa_stats = moa_stats.sort_values("mean_di_resid")

    bottom_10 = moa_stats.head(10)[["moa", "mean_di_resid", "std_di_resid", "n"]].to_dict("records")
    top_10 = moa_stats.tail(10)[["moa", "mean_di_resid", "std_di_resid", "n"]].to_dict("records")

    print(f"\n[{ts()}] LOWEST DI MoAs (most consistent across cell lines):")
    for row in bottom_10:
        print(f"  {row['moa']}: DI_resid = {row['mean_di_resid']:.4f} +/- {row['std_di_resid']:.4f} (n={row['n']})")

    print(f"\n[{ts()}] HIGHEST DI MoAs (most context-dependent):")
    for row in top_10:
        print(f"  {row['moa']}: DI_resid = {row['mean_di_resid']:.4f} +/- {row['std_di_resid']:.4f} (n={row['n']})")

    result = {
        "experiment": "F_moa_stratification",
        "n_drugs": len(moa_sub),
        "n_moa_classes": len(big_moas),
        "min_class_size": MIN_CLASS_SIZE,
        "raw_eta_squared": round(raw_eta, 4),
        "residualized_eta_squared": round(resid_eta, 4),
        "confounds_removed": ["n_cell_lines", "mean_norm"],
        "f_statistic": round(f_stat, 2),
        "f_p_value": f_pval,
        "permutation_n": N_PERMUTATIONS,
        "permutation_null_mean": round(null_etas.mean(), 4),
        "permutation_null_sd": round(null_etas.std(), 4),
        "permutation_null_95th": round(float(np.percentile(null_etas, 95)), 4),
        "permutation_null_99th": round(float(np.percentile(null_etas, 99)), 4),
        "permutation_p": perm_p,
        "alpha": ALPHA,
        "passes": passes,
        "lowest_di_moas": bottom_10,
        "highest_di_moas": top_10,
        "framing": "characterization (DI is mechanism-structured), not prediction (MoA predicts DI)",
        "circularity_note": "MoA and DI may share structural causes (similar targets -> similar cross-context behavior); eta^2 should not be interpreted as causal",
    }

    print(f"\n[{ts()}] {'PASS' if passes else 'FAIL'}")
    save_results(result, RESULTS_PATH)


if __name__ == "__main__":
    main()
