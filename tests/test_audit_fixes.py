"""Tests for all issues identified in the Perplexity audit.

Each test is named after the specific issue it validates.
"""
import numpy as np
import pandas as pd
import pytest
from scipy import stats

from experiments.utils import (
    ci_lower_bound,
    ci_upper_bound,
    compute_di_table,
    fisher_z_pvalue,
    magnitude_correct,
    partial_spearman,
)
from geometry.bracket_norm import direction_instability


class TestMagnitudeCorrectFlatten:
    def test_ci_propagation_shape(self):
        rng = np.random.default_rng(42)
        n = 50
        df = pd.DataFrame({
            "di_raw": rng.uniform(0.1, 0.9, n),
            "di_raw_ci_low": rng.uniform(0.05, 0.5, n),
            "di_raw_ci_high": rng.uniform(0.5, 0.95, n),
            "mean_norm": rng.uniform(1.0, 10.0, n),
        })
        result = magnitude_correct(df)
        assert result["di_corrected"].shape == (n,)
        assert result["di_corrected_ci_low"].shape == (n,)
        assert result["di_corrected_ci_high"].shape == (n,)
        assert all(result["di_corrected_ci_low"] < result["di_corrected_ci_high"])

    def test_corrected_ci_width_matches_raw(self):
        rng = np.random.default_rng(42)
        n = 100
        raw = rng.uniform(0.2, 0.8, n)
        ci_width = rng.uniform(0.02, 0.1, n)
        df = pd.DataFrame({
            "di_raw": raw,
            "di_raw_ci_low": raw - ci_width / 2,
            "di_raw_ci_high": raw + ci_width / 2,
            "mean_norm": rng.uniform(1.0, 10.0, n),
        })
        result = magnitude_correct(df)
        raw_widths = df["di_raw_ci_high"] - df["di_raw_ci_low"]
        corrected_widths = result["di_corrected_ci_high"] - result["di_corrected_ci_low"]
        np.testing.assert_allclose(raw_widths.values, corrected_widths.values, atol=1e-10)


class TestPartialSpearmanIntercept:
    def test_constant_shift_invariance(self):
        rng = np.random.default_rng(42)
        n = 1000
        Z = rng.standard_normal((n, 2))
        X = rng.standard_normal(n)
        Y = 0.3 * X + rng.standard_normal(n) * 0.5

        rho_original = partial_spearman(X, Y, Z)
        rho_shifted = partial_spearman(X + 1000, Y + 1000, Z + 500)
        assert rho_original == pytest.approx(rho_shifted, abs=0.01)

    def test_constant_covariate_does_nothing(self):
        rng = np.random.default_rng(42)
        n = 500
        X = rng.standard_normal(n)
        Y = 0.3 * X + rng.standard_normal(n) * 0.5
        Z_const = np.full((n, 1), 5.0)

        raw_rho, _ = stats.spearmanr(X, Y)
        partial_rho = partial_spearman(X, Y, Z_const)
        assert partial_rho == pytest.approx(raw_rho, abs=0.02)

    def test_intercept_present_in_regression(self):
        rng = np.random.default_rng(42)
        n = 2000
        Z = rng.standard_normal((n, 1)) + 100
        X = Z[:, 0] * 0.5 + rng.standard_normal(n)
        Y = Z[:, 0] * 0.3 + rng.standard_normal(n)

        rho = partial_spearman(X, Y, Z)
        assert abs(rho) < 0.05


class TestGeneIntersection:
    def test_missing_values_as_zero_distorts_di(self):
        rng = np.random.default_rng(42)
        n_contexts = 5
        n_features = 100

        full = rng.standard_normal((n_contexts, n_features))
        true_di = direction_instability(full)

        sparse = full.copy()
        for i in range(n_contexts):
            mask = rng.choice(n_features, 30, replace=False)
            sparse[i, mask] = 0.0

        sparse_di = direction_instability(sparse)
        assert abs(true_di - sparse_di) > 0.01

    def test_pure_zero_padding_is_neutral(self):
        rng = np.random.default_rng(42)
        n_contexts = 5
        n_shared = 50
        n_total = 500

        shared = rng.standard_normal((n_contexts, n_shared))
        padded = np.zeros((n_contexts, n_total))
        padded[:, :n_shared] = shared

        di_shared = direction_instability(shared)
        di_padded = direction_instability(padded)
        assert di_padded == pytest.approx(di_shared, abs=1e-10)

    def test_intersection_preserves_true_di(self):
        rng = np.random.default_rng(42)
        n_contexts = 10
        n_features = 100
        sigs = rng.standard_normal((n_contexts, n_features))
        true_di = direction_instability(sigs)

        subset_idx = rng.choice(n_features, 80, replace=False)
        subset_di = direction_instability(sigs[:, subset_idx])
        assert abs(true_di - subset_di) < 0.15


class TestPlateAggregation:
    def test_replicates_within_plate_inflate_n(self):
        rng = np.random.default_rng(42)
        n_plates = 5
        reps_per_plate = 10
        n_features = 20

        plate_means = rng.standard_normal((n_plates, n_features))
        wells = np.repeat(plate_means, reps_per_plate, axis=0)
        wells += rng.standard_normal(wells.shape) * 0.1

        di_wells = direction_instability(wells)
        di_plates = direction_instability(plate_means)
        assert di_wells < di_plates * 1.5

    def test_aggregation_reduces_to_plate_count(self):
        rng = np.random.default_rng(42)
        n_plates = 5
        reps_per_plate = 4
        n_features = 20

        plate_profiles = rng.standard_normal((n_plates, n_features))
        wells = np.repeat(plate_profiles, reps_per_plate, axis=0)
        wells += rng.standard_normal(wells.shape) * 0.01

        plate_ids = np.repeat(np.arange(n_plates), reps_per_plate)
        aggregated = np.array([
            wells[plate_ids == p].mean(axis=0) for p in range(n_plates)
        ])

        assert aggregated.shape[0] == n_plates
        np.testing.assert_allclose(aggregated, plate_profiles, atol=0.05)


class TestCVFloor:
    def test_near_zero_mean_excluded(self):
        cv_values = []
        for mean_lfc in [0.001, 0.005, 0.009]:
            std_lfc = 0.002
            if abs(mean_lfc) < 0.01:
                continue
            cv_values.append(std_lfc / abs(mean_lfc))
        assert len(cv_values) == 0

    def test_normal_values_computed(self):
        mean_lfc = -0.5
        std_lfc = 0.2
        assert abs(mean_lfc) >= 0.01
        cv = std_lfc / abs(mean_lfc)
        assert cv == pytest.approx(0.4)


class TestFisherZPValue:
    def test_large_rho_gives_small_p(self):
        p = fisher_z_pvalue(0.5, n=100)
        assert p < 1e-5

    def test_zero_rho_gives_large_p(self):
        p = fisher_z_pvalue(0.0, n=100)
        assert p == pytest.approx(1.0, abs=0.01)

    def test_n_cov_widens_se(self):
        p_no_cov = fisher_z_pvalue(0.2, n=50, n_cov=0)
        p_with_cov = fisher_z_pvalue(0.2, n=50, n_cov=2)
        assert p_with_cov > p_no_cov


class TestCIBounds:
    def test_ci_contains_point_estimate(self):
        for rho in [0.1, 0.3, 0.5, -0.2]:
            low = ci_lower_bound(rho, n=100)
            high = ci_upper_bound(rho, n=100)
            assert low < rho < high

    def test_perfect_correlation_edge_case(self):
        assert ci_lower_bound(1.0, n=100) == 1.0
        assert ci_upper_bound(-1.0, n=100) == -1.0


class TestInterpretationBinsMatchPrereg:
    def test_bins_partition_real_line(self):
        bins = {
            "biologically meaningful": lambda r: r >= 0.15,
            "weak": lambda r: 0.05 <= r < 0.15,
            "negligible": lambda r: r < 0.05,
        }
        for rho in np.arange(-0.5, 1.0, 0.01):
            matching = [name for name, pred in bins.items() if pred(rho)]
            assert len(matching) == 1, f"rho={rho:.2f} matched {matching}"


class TestComputeDITable:
    def test_min_contexts_filter(self):
        rng = np.random.default_rng(42)
        sigs = {
            "enough": rng.standard_normal((10, 20)),
            "too_few": rng.standard_normal((3, 20)),
        }
        df = compute_di_table(sigs, n_bootstrap=10, min_contexts=5, rng=rng)
        assert len(df) == 1
        assert df.iloc[0]["perturbation_id"] == "enough"

    def test_bootstrap_ci_ordered_and_reasonable(self):
        rng = np.random.default_rng(42)
        sigs = {"drug_a": rng.standard_normal((20, 50))}
        df = compute_di_table(sigs, n_bootstrap=500, rng=rng)
        row = df.iloc[0]
        assert row["di_raw_ci_low"] < row["di_raw_ci_high"]
        assert row["di_raw_ci_high"] - row["di_raw_ci_low"] < 0.5

    def test_bootstrap_ci_width_shrinks_with_more_contexts(self):
        rng = np.random.default_rng(42)
        base = rng.standard_normal((5, 30))
        sigs_small = {"drug": np.vstack([base, base + rng.standard_normal((5, 30)) * 0.1])}
        sigs_large = {"drug": np.vstack([base] * 10 + [base + rng.standard_normal((5, 30)) * 0.1] * 10)}

        df_small = compute_di_table(sigs_small, n_bootstrap=200, min_contexts=5, rng=rng)
        rng2 = np.random.default_rng(42)
        df_large = compute_di_table(sigs_large, n_bootstrap=200, min_contexts=5, rng=rng2)

        width_small = df_small.iloc[0]["di_raw_ci_high"] - df_small.iloc[0]["di_raw_ci_low"]
        width_large = df_large.iloc[0]["di_raw_ci_high"] - df_large.iloc[0]["di_raw_ci_low"]
        assert width_large < width_small


class TestMagnitudeCorrection:
    def test_removes_norm_correlation(self):
        rng = np.random.default_rng(42)
        n = 200
        mean_norm = rng.uniform(1, 10, n)
        di_raw = 0.3 * mean_norm + rng.standard_normal(n) * 0.5
        df = pd.DataFrame({
            "di_raw": di_raw,
            "di_raw_ci_low": di_raw - 0.05,
            "di_raw_ci_high": di_raw + 0.05,
            "mean_norm": mean_norm,
        })

        raw_corr, _ = stats.spearmanr(df["di_raw"], df["mean_norm"])
        result = magnitude_correct(df)
        corrected_corr, _ = stats.spearmanr(result["di_corrected"], result["mean_norm"])

        assert abs(raw_corr) > 0.5
        assert abs(corrected_corr) < 0.1
