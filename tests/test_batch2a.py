"""Tests for Batch 2A experiments (E2, E3, E5).

Validates statistical methods, magnitude corrections, and decision
criteria match the pre-registration.
"""
import numpy as np
import pandas as pd
import pytest
from scipy import stats
from sklearn.linear_model import LinearRegression

from experiments.utils import (
    ci_lower_bound,
    ci_upper_bound,
    fisher_z_pvalue,
    gini_coefficient,
    partial_spearman,
)


class TestGiniCoefficient:
    def test_uniform_distribution_gini_zero(self):
        values = np.ones(100)
        assert gini_coefficient(values) == pytest.approx(0.0, abs=1e-10)

    def test_maximally_unequal_gini_near_one(self):
        values = np.zeros(1000)
        values[0] = 1.0
        g = gini_coefficient(values)
        assert g > 0.99

    def test_known_value(self):
        values = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        g = gini_coefficient(values)
        expected = 0.2667
        assert g == pytest.approx(expected, abs=0.01)

    def test_empty_returns_zero(self):
        assert gini_coefficient(np.array([])) == 0.0


class TestGiniMagnitudeCorrection:
    def test_correction_removes_magnitude_correlation(self):
        rng = np.random.default_rng(42)
        n = 500
        mean_kill = rng.uniform(0.1, 3.0, n)
        gini_raw = 0.4 * mean_kill + rng.standard_normal(n) * 0.2

        raw_corr, _ = stats.spearmanr(gini_raw, mean_kill)
        assert abs(raw_corr) > 0.5

        X = mean_kill.reshape(-1, 1)
        reg = LinearRegression().fit(X, gini_raw)
        gini_corrected = gini_raw - reg.predict(X) + reg.intercept_

        corrected_corr, _ = stats.spearmanr(gini_corrected, mean_kill)
        assert abs(corrected_corr) < 0.1

    def test_correction_preserves_relative_order_orthogonal_to_magnitude(self):
        rng = np.random.default_rng(42)
        n = 200
        mean_kill = np.ones(n) * 1.5
        true_selectivity = rng.standard_normal(n)
        gini_raw = 0.5 + true_selectivity * 0.1

        X = mean_kill.reshape(-1, 1)
        reg = LinearRegression().fit(X, gini_raw)
        gini_corrected = gini_raw - reg.predict(X) + reg.intercept_

        rho_raw, _ = stats.spearmanr(true_selectivity, gini_raw)
        rho_corr, _ = stats.spearmanr(true_selectivity, gini_corrected)
        assert rho_corr == pytest.approx(rho_raw, abs=0.05)


class TestNoDoubleCounting:
    def test_e3_uses_one_covariate_not_two(self):
        rng = np.random.default_rng(42)
        n = 300
        n_ctx = rng.integers(5, 50, n).astype(float)
        mean_kill = rng.uniform(0.1, 3.0, n)
        gini_raw = 0.3 * mean_kill + rng.standard_normal(n) * 0.2

        X_mag = mean_kill.reshape(-1, 1)
        reg = LinearRegression().fit(X_mag, gini_raw)
        gini_corrected = gini_raw - reg.predict(X_mag) + reg.intercept_

        di = rng.standard_normal(n)

        rho_one_cov = partial_spearman(
            di, gini_corrected, n_ctx.reshape(-1, 1)
        )
        rho_two_cov = partial_spearman(
            di, gini_corrected,
            np.column_stack([n_ctx, mean_kill]),
        )

        assert abs(rho_one_cov - rho_two_cov) < 0.1


class TestBatch2AAlpha:
    def test_alpha_is_0167_not_0125(self):
        alpha = 0.05 / 3
        assert alpha == pytest.approx(0.0167, abs=0.001)

        alpha_batch1 = 0.05 / 4
        assert alpha_batch1 == pytest.approx(0.0125, abs=0.001)
        assert alpha != pytest.approx(alpha_batch1, abs=0.001)


class TestDecisionCriteriaBothConditions:
    def test_ci_gate_blocks_trivial_effects(self):
        n = 7000
        rho = 0.03
        ci_low = ci_lower_bound(rho, n, n_cov=1)
        p = fisher_z_pvalue(rho, n, n_cov=1)

        assert p < 0.0167
        assert ci_low < 0.05
        passes = ci_low > 0.05 and p < 0.0167
        assert not passes

    def test_real_effect_passes_both(self):
        n = 7000
        rho = 0.10
        ci_low = ci_lower_bound(rho, n, n_cov=1)
        p = fisher_z_pvalue(rho, n, n_cov=1)

        assert p < 0.0167
        assert ci_low > 0.05
        passes = ci_low > 0.05 and p < 0.0167
        assert passes

    def test_e5b_negative_criterion(self):
        n = 5000
        rho = -0.10
        ci_high = ci_upper_bound(rho, n, n_cov=2)
        p = fisher_z_pvalue(rho, n, n_cov=2)

        assert ci_high < -0.05
        assert p < 0.0167
        passes = ci_high < -0.05 and p < 0.0167
        assert passes


class TestPartialSpearmanMethodConsistency:
    def test_null_false_positive_rate_below_alpha(self):
        rng = np.random.default_rng(42)
        alpha = 0.0167
        n = 5000
        n_sims = 500
        fp_count = 0

        for _ in range(n_sims):
            X = rng.standard_normal(n)
            Y = rng.standard_normal(n)
            Z = rng.standard_normal((n, 1))

            rho = partial_spearman(X, Y, Z)
            ci_low = ci_lower_bound(rho, n, n_cov=1)
            p = fisher_z_pvalue(rho, n, n_cov=1)

            if ci_low > 0.05 and p < alpha:
                fp_count += 1

        fpr = fp_count / n_sims
        assert fpr < 0.02

    def test_rank_then_residualize_inflates_with_correlated_covariate(self):
        rng = np.random.default_rng(99)
        n = 5000
        correct_rhos = []
        wrong_rhos = []
        n_sims = 200

        for _ in range(n_sims):
            Z = rng.standard_normal((n, 2))
            X = Z @ rng.standard_normal(2) * 0.3 + rng.standard_normal(n)
            Y = Z @ rng.standard_normal(2) * 0.3 + rng.standard_normal(n)

            correct_rhos.append(partial_spearman(X, Y, Z))

            X_r = stats.rankdata(X)
            Y_r = stats.rankdata(Y)
            Z_r = np.column_stack([stats.rankdata(Z[:, i]) for i in range(2)])
            X_res = X_r - Z_r @ np.linalg.lstsq(Z_r, X_r, rcond=None)[0]
            Y_res = Y_r - Z_r @ np.linalg.lstsq(Z_r, Y_r, rcond=None)[0]
            wrong_rhos.append(np.corrcoef(X_res, Y_res)[0, 1])

        assert abs(np.mean(correct_rhos)) < 0.02
        assert np.mean(wrong_rhos) > 0.05


class TestE5bExpectedSign:
    def test_broad_compensation_lowers_di(self):
        rng = np.random.default_rng(42)
        n = 1000
        paralog_breadth = rng.uniform(0, 1, n)
        di = -0.15 * paralog_breadth + rng.standard_normal(n) * 0.5
        Z = rng.standard_normal((n, 2))

        rho = partial_spearman(di, paralog_breadth, Z)
        assert rho < 0
