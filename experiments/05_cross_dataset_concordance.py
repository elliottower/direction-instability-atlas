"""Experiment A: Cross-dataset DI concordance (A1 + A2).

A1: LINCS ∩ Tahoe (n=103, same modality)
A2: LINCS ∩ JUMP-CP compounds (n=330, cross-modality)

For each pair: raw Spearman rho, partial Spearman rho controlling
for n_contexts in both datasets (2 covariates for A1; 1 covariate
for A2 since JUMP-CP compound context counts are unavailable in
the pre-computed npz), 95% CI via Fisher z.

Decision criterion: CI lower bound > 0.05, p < 0.0125.
Partial Spearman method: Pearson-residualize then Spearman.

Inputs:
  - Pre-computed LINCS DI: drug-perturbation-geometry/zenodo_v1/drug_instability_8949.csv
  - Tahoe DI: results/tahoe_di.csv (from script 01)
  - Pre-computed JUMP-CP compound DI: bracket-norm-validity/results/08_jump_cp/
  - JUMP-CP compound metadata: bracket-norm-validity/data/jump_cp/compound_metadata.csv.gz

Output: results/concordance.json
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
    magnitude_correct,
    partial_spearman,
    save_results,
)

LINCS_PATH = Path("../drug-perturbation-geometry/zenodo_v1/drug_instability_8949.csv")
TAHOE_PATH = Path("results/tahoe_di.csv")
JUMP_DI_PATH = Path("../bracket-norm-validity/results/08_jump_cp/jump_cp_instabilities.npz")
JUMP_META_PATH = Path("../bracket-norm-validity/data/jump_cp/compound_metadata.csv.gz")
RESULTS_PATH = Path("results/concordance.json")

ALPHA = 0.0125
CI_FLOOR = 0.05


def load_lincs_di() -> pd.DataFrame:
    df = pd.read_csv(LINCS_PATH)
    df = df.rename(columns={
        "drug_name": "drug",
        "direction_instability": "di_raw",
        "n_cell_lines": "n_contexts",
    })
    df = df[["drug", "di_raw", "n_contexts", "mean_norm"]]
    df["di_raw_ci_low"] = df["di_raw"]
    df["di_raw_ci_high"] = df["di_raw"]
    df = magnitude_correct(df)
    return df


def load_tahoe_di() -> pd.DataFrame:
    df = pd.read_csv(TAHOE_PATH)
    return df[["drug", "di_corrected", "di_raw", "n_contexts"]]


def load_jump_compound_di() -> pd.DataFrame:
    data = np.load(JUMP_DI_PATH, allow_pickle=True)
    meta = pd.read_csv(JUMP_META_PATH)
    jcp_to_inchikey = dict(zip(meta["Metadata_JCP2022"], meta["Metadata_InChIKey"]))

    df = pd.DataFrame({
        "compound_id": data["compounds"],
        "di_raw": data["raw"],
        "di_corrected": data["corrected"],
    })
    df["inchikey"] = df["compound_id"].map(jcp_to_inchikey)
    return df


def test_concordance(
    name: str,
    di_a: np.ndarray,
    di_b: np.ndarray,
    covariates: np.ndarray,
    n_cov: int,
) -> dict:
    n = len(di_a)

    raw_rho, raw_p = stats.spearmanr(di_a, di_b)

    partial_rho = partial_spearman(di_a, di_b, covariates)

    ci_low = ci_lower_bound(partial_rho, n, n_cov=n_cov)
    ci_high = ci_upper_bound(partial_rho, n, n_cov=n_cov)
    p_val = fisher_z_pvalue(partial_rho, n, n_cov=n_cov)

    passes = ci_low > CI_FLOOR and p_val < ALPHA

    return {
        "test": name,
        "n": n,
        "n_covariates": n_cov,
        "raw_spearman_rho": round(raw_rho, 4),
        "raw_p_value": raw_p,
        "partial_spearman_rho": round(partial_rho, 4),
        "partial_ci_low": round(ci_low, 4),
        "partial_ci_high": round(ci_high, 4),
        "partial_p_value": p_val,
        "criterion": f"CI lower > {CI_FLOOR}",
        "alpha": ALPHA,
        "passes": passes,
    }


def main():
    ts = lambda: datetime.now().strftime("%H:%M:%S")
    print(f"[{ts()}] Experiment A: Cross-dataset DI concordance")

    lincs = load_lincs_di()
    tahoe = load_tahoe_di()
    print(f"[{ts()}] LINCS: {len(lincs)} drugs, Tahoe: {len(tahoe)} drugs")

    # A1: LINCS ∩ Tahoe (matched by drug name, InChIKey-verified n=103)
    merged_a1 = lincs.merge(tahoe, on="drug", suffixes=("_lincs", "_tahoe"))
    assert len(merged_a1) >= 100, (
        f"Expected ~103 LINCS-Tahoe overlap, got {len(merged_a1)}"
    )
    print(f"[{ts()}] A1 overlap: {len(merged_a1)} drugs")

    Z_a1 = np.column_stack([
        merged_a1["n_contexts_lincs"].values,
        merged_a1["n_contexts_tahoe"].values,
    ]).astype(np.float64)

    result_a1 = test_concordance(
        "A1_LINCS_Tahoe",
        merged_a1["di_corrected_lincs"].values.astype(np.float64),
        merged_a1["di_corrected_tahoe"].values.astype(np.float64),
        Z_a1,
        n_cov=2,
    )
    print(f"[{ts()}] A1: partial rho = {result_a1['partial_spearman_rho']:.4f} "
          f"[{result_a1['partial_ci_low']:.4f}, {result_a1['partial_ci_high']:.4f}] "
          f"{'PASS' if result_a1['passes'] else 'FAIL'}")

    # A2: LINCS ∩ JUMP-CP compounds (matched by InChIKey connectivity)
    jump = load_jump_compound_di()
    print(f"[{ts()}] JUMP-CP: {len(jump)} compounds")

    lincs_pert_info_path = Path(
        "../drug-perturbation-geometry/data/GSE92742_Broad_LINCS_pert_info.txt.gz"
    )
    if lincs_pert_info_path.exists():
        lincs_meta = pd.read_csv(lincs_pert_info_path, sep="\t")
        lincs_meta = lincs_meta[lincs_meta["pert_type"] == "trt_cp"]
        lincs_meta["inchikey_prefix"] = lincs_meta["inchi_key"].str[:14]

        lincs_with_ik = lincs.merge(
            lincs_meta[["pert_iname", "inchikey_prefix"]].drop_duplicates(),
            left_on="drug", right_on="pert_iname", how="inner",
        )

        jump["inchikey_prefix"] = jump["inchikey"].str[:14]

        merged_a2 = lincs_with_ik.merge(
            jump, on="inchikey_prefix", suffixes=("_lincs", "_jump"),
        )
        print(f"[{ts()}] A2 overlap (InChIKey): {len(merged_a2)} compounds")

        merged_a2 = merged_a2.drop_duplicates(subset="drug")

        if len(merged_a2) > 0:
            Z_a2 = merged_a2["n_contexts"].values.astype(np.float64).reshape(-1, 1)

            di_lincs = merged_a2["di_corrected_lincs"].values.astype(np.float64)
            di_jump = merged_a2["di_corrected_jump"].values.astype(np.float64)

            result_a2 = test_concordance(
                "A2_LINCS_JUMP",
                di_lincs,
                di_jump,
                Z_a2,
                n_cov=1,
            )
            result_a2["note"] = "1 covariate (LINCS n_contexts only); JUMP-CP context counts unavailable in pre-computed npz"
        else:
            result_a2 = {"test": "A2_LINCS_JUMP", "n": 0, "passes": False,
                         "note": "No InChIKey overlap found"}
    else:
        result_a2 = {"test": "A2_LINCS_JUMP", "note": "LINCS pert_info not found"}

    print(f"[{ts()}] A2: {result_a2}")

    results = {
        "experiment": "A_concordance",
        "alpha": ALPHA,
        "ci_floor": CI_FLOOR,
        "method": "Pearson-residualize then Spearman",
        "A1": result_a1,
        "A2": result_a2,
    }
    save_results(results, RESULTS_PATH)


if __name__ == "__main__":
    main()
