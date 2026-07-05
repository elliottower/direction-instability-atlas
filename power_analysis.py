"""Power analysis for pre-registration thresholds.

Uses CORRECT partial Spearman rho method: Pearson-residualize raw values,
then compute Spearman on the residuals. The "rank then residualize" method
is biased (gives ~0.25 when true = 0).

Criterion: 95% CI lower bound > floor (not point estimate > threshold).
CI via Fisher z-transform on the partial correlation.
"""
import json
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy import stats


def partial_spearman(X: np.ndarray, Y: np.ndarray, Z: np.ndarray) -> float:
    """Correct partial Spearman rho: Pearson-residualize, then Spearman."""
    X_resid = X - Z @ np.linalg.lstsq(Z, X, rcond=None)[0]
    Y_resid = Y - Z @ np.linalg.lstsq(Z, Y, rcond=None)[0]
    return stats.spearmanr(X_resid, Y_resid)[0]


def ci_lower_bound(rho: float, n: int, n_cov: int) -> float:
    """Fisher z 95% CI lower bound for partial correlation."""
    if abs(rho) >= 1.0:
        return rho
    z = np.arctanh(rho)
    se = 1.0 / np.sqrt(n - n_cov - 3)
    return float(np.tanh(z - 1.96 * se))


def simulate_power_ci(
    n: int,
    true_rho: float,
    ci_floor: float,
    n_covariates: int = 2,
    n_sims: int = 5000,
    alpha: float = 0.0125,
    rng: np.random.Generator | None = None,
) -> dict:
    """Power = P(CI lower bound > floor AND p < alpha)."""
    if rng is None:
        rng = np.random.default_rng()

    passes = 0
    observed_rhos = []

    for _ in range(n_sims):
        Z = rng.standard_normal((n, n_covariates))
        noise_x = rng.standard_normal(n)
        noise_y = rng.standard_normal(n)
        coef_x = rng.standard_normal(n_covariates) * 0.3
        coef_y = rng.standard_normal(n_covariates) * 0.3
        X = Z @ coef_x + noise_x
        Y = Z @ coef_y + true_rho * noise_x + np.sqrt(1 - true_rho**2) * noise_y

        rho = partial_spearman(X, Y, Z)
        observed_rhos.append(rho)

        ci_low = ci_lower_bound(rho, n, n_covariates)
        z = np.arctanh(rho) if abs(rho) < 1.0 else 10.0
        se = 1.0 / np.sqrt(n - n_covariates - 3)
        p = 2 * (1 - stats.norm.cdf(abs(z) / se))

        if ci_low > ci_floor and p < alpha:
            passes += 1

    observed_rhos = np.array(observed_rhos)
    return {
        "n": n,
        "true_rho": true_rho,
        "ci_floor": ci_floor,
        "alpha": alpha,
        "n_covariates": n_covariates,
        "n_sims": n_sims,
        "power": passes / n_sims,
        "mean_obs": float(np.mean(observed_rhos)),
        "std_obs": float(np.std(observed_rhos)),
    }


def simulate_power_d(
    n: int,
    true_rho: float,
    ci_floor: float,
    n_sims: int = 5000,
    alpha: float = 0.0125,
    rng: np.random.Generator | None = None,
) -> dict:
    """Power for Exp D: CI lower bound > floor, no covariates."""
    if rng is None:
        rng = np.random.default_rng()

    passes = 0
    observed = []

    for _ in range(n_sims):
        noise_x = rng.standard_normal(n)
        noise_y = rng.standard_normal(n)
        X = noise_x
        Y = true_rho * noise_x + np.sqrt(1 - true_rho**2) * noise_y

        rho, p = stats.spearmanr(X, Y)
        observed.append(rho)

        cl = ci_lower_bound(rho, n, 0)
        if cl > ci_floor and p < alpha:
            passes += 1

    observed = np.array(observed)
    return {
        "n": n,
        "true_rho": true_rho,
        "ci_floor": ci_floor,
        "alpha": alpha,
        "n_sims": n_sims,
        "power": passes / n_sims,
        "mean_obs": float(np.mean(observed)),
    }


def main():
    ts = lambda: datetime.now().strftime("%H:%M:%S")
    print(f"[{ts()}] Power analysis — CORRECTED partial Spearman method")
    print(f"[{ts()}] Method: Pearson-residualize then Spearman on residuals")
    print(f"[{ts()}] Criterion: 95% CI lower bound > floor, alpha = 0.0125 (4 tests)")
    print("=" * 70)

    rng = np.random.default_rng(2026)
    all_results = {}

    # === Sanity check: null should give ~0 ===
    print(f"\n[{ts()}] Sanity check: true_rho=0, n=1000")
    r_null = simulate_power_ci(n=1000, true_rho=0.0, ci_floor=0.10, n_sims=2000, rng=rng)
    print(f"  mean observed partial rho = {r_null['mean_obs']:.4f} (should be ~0)")
    print(f"  FPR (CI lower > 0.10) = {r_null['power']:.4f} (should be ~0)")

    # === A1: LINCS ∩ Tahoe, n=103, 2 covariates ===
    print(f"\n[{ts()}] === Experiment A1: LINCS ∩ Tahoe (n=103) ===")
    print("  CI-lower-bound > 0.05, alpha = 0.0125")
    a1 = []
    for true_rho in [0.20, 0.25, 0.30, 0.35, 0.40]:
        r = simulate_power_ci(n=103, true_rho=true_rho, ci_floor=0.05, rng=rng)
        a1.append(r)
        print(f"  true_rho={true_rho:.2f}: power={r['power']:.3f}, mean_obs={r['mean_obs']:.3f}")
    all_results["A1"] = a1

    # === A2: LINCS ∩ JUMP-CP, n=330, 2 covariates ===
    print(f"\n[{ts()}] === Experiment A2: LINCS ∩ JUMP-CP (n=330) ===")
    print("  CI-lower-bound > 0.05, alpha = 0.0125")
    a2 = []
    for true_rho in [0.10, 0.15, 0.20, 0.25, 0.30]:
        r = simulate_power_ci(n=330, true_rho=true_rho, ci_floor=0.05, rng=rng)
        a2.append(r)
        print(f"  true_rho={true_rho:.2f}: power={r['power']:.3f}, mean_obs={r['mean_obs']:.3f}")
    all_results["A2"] = a2

    # === C: CRISPR essentiality, n=7948, 1 covariate ===
    print(f"\n[{ts()}] === Experiment C: CRISPR DI vs essentiality (n=7948) ===")
    print("  CI-lower-bound > 0.05, alpha = 0.0125")
    c_results = []
    for true_rho in [0.05, 0.08, 0.10, 0.15, 0.20]:
        r = simulate_power_ci(n=7948, true_rho=true_rho, ci_floor=0.05,
                              n_covariates=1, n_sims=3000, rng=rng)
        c_results.append(r)
        print(f"  true_rho={true_rho:.2f}: power={r['power']:.3f}, mean_obs={r['mean_obs']:.3f}")
    all_results["C"] = c_results

    # === D: CRISPR vs ORF, n=5220, no covariates, CI-based ===
    print(f"\n[{ts()}] === Experiment D: CRISPR vs ORF DI (n=5220) ===")
    print("  CI lower > 0.05, alpha = 0.0125 (no covariates, CI-based)")
    d_results = []
    for true_rho in [0.05, 0.08, 0.10, 0.15, 0.20]:
        r = simulate_power_d(n=5220, true_rho=true_rho, ci_floor=0.05,
                             n_sims=3000, rng=rng)
        d_results.append(r)
        print(f"  true_rho={true_rho:.2f}: power={r['power']:.3f}, mean_obs={r['mean_obs']:.3f}")
    all_results["D"] = d_results

    # === FPR check for each experiment ===
    print(f"\n[{ts()}] === False positive rates (true_rho = 0) ===")
    for label, n, ci_floor, n_cov in [
        ("A1", 103, 0.05, 2),
        ("A2", 330, 0.05, 2),
        ("C", 7948, 0.05, 1),
    ]:
        r = simulate_power_ci(n=n, true_rho=0.0, ci_floor=ci_floor,
                              n_covariates=n_cov, n_sims=3000, rng=rng)
        print(f"  {label} (n={n}): FPR = {r['power']:.4f}, mean_obs = {r['mean_obs']:.4f}")

    r_d = simulate_power_d(n=5220, true_rho=0.0, ci_floor=0.05,
                           n_sims=3000, rng=rng)
    print(f"  D  (n=5220): FPR = {r_d['power']:.4f}, mean_obs = {r_d['mean_obs']:.4f}")

    # Save
    output_path = Path("data/power_analysis_results.json")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump({
            "generated": datetime.now().isoformat(),
            "method": "Pearson-residualize then Spearman (corrected)",
            "alpha": 0.0125,
            "n_primary_tests": 4,
            "results": all_results,
        }, f, indent=2)

    print(f"\n[{ts()}] Saved to {output_path}")
    print(f"[{ts()}] Done.")


if __name__ == "__main__":
    main()
