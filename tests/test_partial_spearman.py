"""Tests for the corrected partial Spearman rho computation.

Validates that:
1. The corrected method (Pearson-residualize then Spearman) gives unbiased estimates
2. The buggy method (rank then residualize) is demonstrably biased
3. CI lower bound computation is correct
4. FPR is controlled under the null
"""
import numpy as np
import pytest
from scipy import stats

from power_analysis import ci_lower_bound, partial_spearman


class TestPartialSpearman:
    def test_null_gives_zero(self):
        rng = np.random.default_rng(42)
        n = 2000
        observed = []
        for _ in range(500):
            Z = rng.standard_normal((n, 2))
            X = Z @ rng.standard_normal(2) * 0.3 + rng.standard_normal(n)
            Y = Z @ rng.standard_normal(2) * 0.3 + rng.standard_normal(n)
            observed.append(partial_spearman(X, Y, Z))
        mean_rho = np.mean(observed)
        assert mean_rho == pytest.approx(0.0, abs=0.02)

    def test_known_effect_recovered(self):
        rng = np.random.default_rng(123)
        n = 5000
        true_rho = 0.30
        observed = []
        for _ in range(200):
            Z = rng.standard_normal((n, 2))
            noise_x = rng.standard_normal(n)
            X = Z @ rng.standard_normal(2) * 0.3 + noise_x
            Y = (Z @ rng.standard_normal(2) * 0.3
                 + true_rho * noise_x
                 + np.sqrt(1 - true_rho**2) * rng.standard_normal(n))
            observed.append(partial_spearman(X, Y, Z))
        mean_rho = np.mean(observed)
        assert mean_rho == pytest.approx(true_rho, abs=0.02)

    def test_rank_then_residualize_is_biased(self):
        """Regression test: the old method gives inflated partial rho under null."""
        rng = np.random.default_rng(99)
        n = 5000
        buggy_results = []
        correct_results = []
        for _ in range(200):
            Z = rng.standard_normal((n, 2))
            X = Z @ rng.standard_normal(2) * 0.3 + rng.standard_normal(n)
            Y = Z @ rng.standard_normal(2) * 0.3 + rng.standard_normal(n)

            correct_results.append(partial_spearman(X, Y, Z))

            X_r = stats.rankdata(X)
            Y_r = stats.rankdata(Y)
            Z_r = np.column_stack([stats.rankdata(Z[:, i]) for i in range(2)])
            X_res = X_r - Z_r @ np.linalg.lstsq(Z_r, X_r, rcond=None)[0]
            Y_res = Y_r - Z_r @ np.linalg.lstsq(Z_r, Y_r, rcond=None)[0]
            buggy_results.append(np.corrcoef(X_res, Y_res)[0, 1])

        assert abs(np.mean(correct_results)) < 0.02
        assert np.mean(buggy_results) > 0.05

    def test_single_covariate(self):
        rng = np.random.default_rng(77)
        n = 3000
        true_rho = 0.20
        observed = []
        for _ in range(200):
            Z = rng.standard_normal((n, 1))
            noise_x = rng.standard_normal(n)
            X = Z[:, 0] * 0.4 + noise_x
            Y = (Z[:, 0] * 0.3
                 + true_rho * noise_x
                 + np.sqrt(1 - true_rho**2) * rng.standard_normal(n))
            observed.append(partial_spearman(X, Y, Z))
        assert np.mean(observed) == pytest.approx(true_rho, abs=0.02)


class TestCILowerBound:
    def test_positive_rho_gives_lower_ci(self):
        ci = ci_lower_bound(0.30, n=100, n_cov=2)
        assert ci < 0.30
        assert ci > 0.0

    def test_larger_n_gives_tighter_ci(self):
        ci_small = ci_lower_bound(0.30, n=100, n_cov=2)
        ci_large = ci_lower_bound(0.30, n=1000, n_cov=2)
        assert ci_large > ci_small

    def test_coverage(self):
        """95% CI should cover the true value ~95% of the time."""
        rng = np.random.default_rng(55)
        n = 200
        true_rho = 0.25
        covers = 0
        n_sims = 1000
        for _ in range(n_sims):
            Z = rng.standard_normal((n, 2))
            noise_x = rng.standard_normal(n)
            X = Z @ rng.standard_normal(2) * 0.3 + noise_x
            Y = (Z @ rng.standard_normal(2) * 0.3
                 + true_rho * noise_x
                 + np.sqrt(1 - true_rho**2) * rng.standard_normal(n))
            rho = partial_spearman(X, Y, Z)
            ci_low = ci_lower_bound(rho, n, 2)
            z = np.arctanh(rho) if abs(rho) < 1 else 10
            se = 1.0 / np.sqrt(n - 2 - 3)
            ci_high = float(np.tanh(z + 1.96 * se))
            if ci_low <= true_rho <= ci_high:
                covers += 1
        coverage = covers / n_sims
        assert coverage == pytest.approx(0.95, abs=0.03)


class TestFPRControl:
    def test_fpr_under_alpha(self):
        """Under null, fraction of CI lower > 0.05 AND p < 0.0125 should be < alpha."""
        rng = np.random.default_rng(42)
        n = 300
        alpha = 0.0125
        ci_floor = 0.05
        passes = 0
        n_sims = 2000
        for _ in range(n_sims):
            Z = rng.standard_normal((n, 2))
            X = Z @ rng.standard_normal(2) * 0.3 + rng.standard_normal(n)
            Y = Z @ rng.standard_normal(2) * 0.3 + rng.standard_normal(n)
            rho = partial_spearman(X, Y, Z)
            ci_low = ci_lower_bound(rho, n, 2)
            z = np.arctanh(rho) if abs(rho) < 1 else 10
            se = 1.0 / np.sqrt(n - 2 - 3)
            p = 2 * (1 - stats.norm.cdf(abs(z) / se))
            if ci_low > ci_floor and p < alpha:
                passes += 1
        fpr = passes / n_sims
        assert fpr < alpha
