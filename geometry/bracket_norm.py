"""Bracket norm computation and validity-conditioned variants.

Raw bracket norm measures how much a perturbation response changes DIRECTION
as context varies. This module implements the raw form plus validity-filtered
variants that condition on toxicity, phenotype alignment, and transport stability.

All functions operate on signature matrices (n_contexts, n_features).
"""
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.stats import spearmanr


def direction_instability(signatures: np.ndarray) -> float:
    """Raw bracket norm analogue: 1 - mean(pairwise cosine of unit signatures).

    Args:
        signatures: (n_contexts, n_features) perturbation response vectors.

    Returns:
        Scalar in [0, 2]. Low = stable direction. High = context-dependent.
    """
    norms = np.linalg.norm(signatures, axis=1, keepdims=True)
    norms = np.maximum(norms, 1e-10)
    unit = signatures / norms
    cosines = unit @ unit.T
    n = cosines.shape[0]
    triu_idx = np.triu_indices(n, k=1)
    return 1.0 - float(cosines[triu_idx].mean())


def magnitude_cv(signatures: np.ndarray) -> float:
    """Coefficient of variation of signature magnitudes across contexts."""
    norms = np.linalg.norm(signatures, axis=1)
    if norms.mean() < 1e-10:
        return 0.0
    return float(norms.std() / norms.mean())


def phenotype_projected_bracket(
    signatures: np.ndarray,
    phenotype_direction: np.ndarray,
) -> float:
    """Bracket norm projected onto a phenotype target direction.

    Asks: does the direction instability push TOWARD the desired phenotype,
    or is it orthogonal/opposing noise?

    Args:
        signatures: (n_contexts, n_features) perturbation response vectors.
        phenotype_direction: (n_features,) unit vector of desired phenotype change.

    Returns:
        Mean absolute projection of pairwise difference vectors onto phenotype.
        High = instability is phenotype-aligned (potentially informative).
        Low = instability is orthogonal to phenotype (noise or off-target).
    """
    phenotype_direction = phenotype_direction / (np.linalg.norm(phenotype_direction) + 1e-10)
    n = signatures.shape[0]
    projections = []
    for i in range(n):
        for j in range(i + 1, n):
            diff = signatures[i] - signatures[j]
            proj = abs(float(diff @ phenotype_direction))
            projections.append(proj)
    return float(np.mean(projections)) if projections else 0.0


def toxicity_corrected_bracket(
    signatures: np.ndarray,
    stress_genes: np.ndarray,
) -> float:
    """Direction instability after removing stress/toxicity gene axes.

    Args:
        signatures: (n_contexts, n_features) perturbation response vectors.
        stress_genes: (n_stress,) indices of stress/apoptosis genes to remove.

    Returns:
        Direction instability computed on the residual (non-stress) genes.
    """
    mask = np.ones(signatures.shape[1], dtype=bool)
    mask[stress_genes] = False
    return direction_instability(signatures[:, mask])


def transport_stable_bracket(
    signatures: np.ndarray,
    frechet_penalty: float = 1.0,
) -> float:
    """Bracket norm penalized by cross-context Frechet variance.

    A high raw bracket that is UNSTABLE across contexts should not support
    a global mechanism claim. This penalizes instability that doesn't replicate.

    Args:
        signatures: (n_contexts, n_features) perturbation response vectors.
        frechet_penalty: weight on variance penalty.

    Returns:
        raw_bracket - frechet_penalty * frechet_variance.
        Positive = genuine cross-context deformation.
        Near-zero or negative = unstable/unreliable bracket.
    """
    raw = direction_instability(signatures)

    norms = np.linalg.norm(signatures, axis=1, keepdims=True)
    norms = np.maximum(norms, 1e-10)
    unit = signatures / norms
    frechet_mean = unit.mean(axis=0)
    frechet_mean = frechet_mean / (np.linalg.norm(frechet_mean) + 1e-10)
    deviations = np.arccos(np.clip(unit @ frechet_mean, -1, 1))
    frechet_var = float(np.mean(deviations**2))

    return raw - frechet_penalty * frechet_var


def localization_score(
    signatures: np.ndarray,
    region_mask: np.ndarray,
) -> float:
    """Whether bracket norm is concentrated in specific gene regions.

    High localization = bracket is driven by a specific pathway (good).
    Low localization = bracket is diffuse/global (suspicious for toxicity).

    Args:
        signatures: (n_contexts, n_features) perturbation response vectors.
        region_mask: (n_features,) boolean mask for genes in a specific pathway.

    Returns:
        Ratio of within-region instability to global instability.
        >1 means instability is concentrated in the region.
    """
    global_bracket = direction_instability(signatures)
    region_bracket = direction_instability(signatures[:, region_mask])
    if global_bracket < 1e-10:
        return 0.0
    return region_bracket / global_bracket


def validity_ladder_score(
    raw_bracket: float,
    has_process_view: bool,
    has_vector_field: bool,
    is_localized: bool,
    phenotype_aligned: bool,
    transports: bool,
    predicts_holdout: bool,
) -> dict:
    """Evaluate bracket norm against the seven-rung validity ladder.

    Returns a dict with per-rung pass/fail and an overall validity level.
    """
    rungs = {
        "1_process_view_declared": has_process_view,
        "2_reliable_vector_field": has_vector_field,
        "3_perturbation_estimable": raw_bracket > 0.01,
        "4_localized_to_biology": is_localized,
        "5_phenotype_aligned": phenotype_aligned,
        "6_transports_across_contexts": transports,
        "7_predicts_holdout": predicts_holdout,
    }
    highest_passed = 0
    for i, (name, passed) in enumerate(rungs.items(), 1):
        if passed:
            highest_passed = i
        else:
            break

    return {
        "rungs": rungs,
        "highest_level": highest_passed,
        "valid_for_mechanism_claim": highest_passed >= 5,
        "raw_bracket": raw_bracket,
    }
