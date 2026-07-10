"""Shared utilities for DI atlas experiments.

All experiments use these functions for DI computation, magnitude
correction, bootstrap CIs, and statistical testing.
"""
import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LinearRegression

from geometry.bracket_norm import direction_instability, magnitude_cv


def compute_di_table(
    signatures_by_id: dict[str, np.ndarray],
    n_bootstrap: int = 1000,
    min_contexts: int = 5,
    rng: np.random.Generator | None = None,
) -> pd.DataFrame:
    """Compute DI + bootstrap CIs + magnitude CV for each perturbation.

    Args:
        signatures_by_id: {perturbation_id: (n_contexts, n_features) array}
        n_bootstrap: number of bootstrap resamples for CIs
        min_contexts: minimum contexts required (skip otherwise)
        rng: random number generator for bootstrap

    Returns:
        DataFrame with columns: perturbation_id, di_raw, di_raw_ci_low,
        di_raw_ci_high, magnitude_cv, mean_norm, n_contexts
    """
    if rng is None:
        rng = np.random.default_rng(42)

    rows = []
    for pid, sigs in signatures_by_id.items():
        n_ctx = sigs.shape[0]
        if n_ctx < min_contexts:
            continue

        di = direction_instability(sigs)
        mcv = magnitude_cv(sigs)
        mean_norm = float(np.linalg.norm(sigs, axis=1).mean())

        boot_dis = []
        for _ in range(n_bootstrap):
            idx = rng.choice(n_ctx, size=n_ctx, replace=True)
            boot_dis.append(direction_instability(sigs[idx]))
        boot_dis = np.array(boot_dis)

        rows.append({
            "perturbation_id": pid,
            "di_raw": di,
            "di_raw_ci_low": float(np.percentile(boot_dis, 2.5)),
            "di_raw_ci_high": float(np.percentile(boot_dis, 97.5)),
            "magnitude_cv": mcv,
            "mean_norm": mean_norm,
            "n_contexts": n_ctx,
        })

    return pd.DataFrame(rows)


def magnitude_correct(df: pd.DataFrame) -> pd.DataFrame:
    """Add di_corrected (+ CIs) by residualizing DI against mean_norm."""
    if len(df) == 0:
        df = df.copy()
        for col in ["di_corrected", "di_corrected_ci_low", "di_corrected_ci_high"]:
            df[col] = pd.Series(dtype=float)
        return df

    X = df["mean_norm"].values.reshape(-1, 1)
    y = df["di_raw"].values

    reg = LinearRegression().fit(X, y)
    residuals = y - reg.predict(X)
    corrected = residuals + reg.intercept_

    df = df.copy()
    df["di_corrected"] = corrected

    y_low = df["di_raw_ci_low"].values
    y_high = df["di_raw_ci_high"].values
    pred = reg.predict(X).flatten()
    df["di_corrected_ci_low"] = (y_low - pred) + reg.intercept_
    df["di_corrected_ci_high"] = (y_high - pred) + reg.intercept_

    return df


def partial_spearman(X: np.ndarray, Y: np.ndarray, Z: np.ndarray) -> float:
    """Correct partial Spearman rho: Pearson-residualize, then Spearman."""
    Z_int = np.column_stack([np.ones(len(Z)), Z])
    X_resid = X - Z_int @ np.linalg.lstsq(Z_int, X, rcond=None)[0]
    Y_resid = Y - Z_int @ np.linalg.lstsq(Z_int, Y, rcond=None)[0]
    return stats.spearmanr(X_resid, Y_resid)[0]


def ci_lower_bound(rho: float, n: int, n_cov: int = 0) -> float:
    """Fisher z 95% CI lower bound for (partial) correlation."""
    if abs(rho) >= 1.0:
        return rho
    z = np.arctanh(rho)
    se = 1.0 / np.sqrt(n - n_cov - 3)
    return float(np.tanh(z - 1.96 * se))


def ci_upper_bound(rho: float, n: int, n_cov: int = 0) -> float:
    """Fisher z 95% CI upper bound for (partial) correlation."""
    if abs(rho) >= 1.0:
        return rho
    z = np.arctanh(rho)
    se = 1.0 / np.sqrt(n - n_cov - 3)
    return float(np.tanh(z + 1.96 * se))


def fisher_z_pvalue(rho: float, n: int, n_cov: int = 0) -> float:
    """Two-sided p-value via Fisher z-transform."""
    if abs(rho) >= 1.0:
        return 0.0
    z = np.arctanh(rho)
    se = 1.0 / np.sqrt(n - n_cov - 3)
    return float(2 * (1 - stats.norm.cdf(abs(z) / se)))


def gini_coefficient(values: np.ndarray) -> float:
    """Gini coefficient of an array of non-negative values."""
    values = np.sort(np.abs(values))
    n = len(values)
    if n == 0 or values.sum() == 0:
        return 0.0
    index = np.arange(1, n + 1)
    return float((2 * np.sum(index * values) - (n + 1) * np.sum(values)) / (n * np.sum(values)))


class _NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.bool_):
            return bool(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)


def save_results(results: dict, path: Path) -> None:
    """Save results dict as JSON with timestamp."""
    results["generated"] = datetime.now().isoformat()
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(results, f, indent=2, cls=_NumpyEncoder)
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Saved to {path}")
