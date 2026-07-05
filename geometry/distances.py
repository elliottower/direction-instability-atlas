"""Distance metrics and geometric utilities."""
import numpy as np
from scipy.linalg import svd
from scipy.spatial.distance import pdist, squareform


def grassmannian_distance(U: np.ndarray, V: np.ndarray) -> float:
    """Geodesic distance between two subspaces on the Grassmannian.

    Args:
        U: (d, k) orthonormal basis for subspace 1.
        V: (d, k) orthonormal basis for subspace 2.

    Returns:
        Geodesic distance = sqrt(sum of squared principal angles).
    """
    M = U.T @ V
    svals = np.clip(svd(M, compute_uv=False), 0, 1)
    angles = np.arccos(svals)
    return float(np.sqrt(np.sum(angles**2)))


def subspace_overlap(U: np.ndarray, V: np.ndarray) -> float:
    """Mean squared cosine of principal angles (0=orthogonal, 1=identical)."""
    M = U.T @ V
    svals = svd(M, compute_uv=False)
    return float(np.mean(svals**2))


def pca_subspace(X: np.ndarray, k: int) -> np.ndarray:
    """Extract top-k PCA subspace from data matrix.

    Args:
        X: (n_samples, n_features) centered data.
        k: number of components.

    Returns:
        (n_features, k) orthonormal basis.
    """
    X_centered = X - X.mean(axis=0, keepdims=True)
    _, _, Vt = svd(X_centered, full_matrices=False)
    return Vt[:k].T


def cosine_similarity_matrix(signatures: np.ndarray) -> np.ndarray:
    """Pairwise cosine similarity between rows of a matrix."""
    norms = np.linalg.norm(signatures, axis=1, keepdims=True)
    norms = np.maximum(norms, 1e-10)
    unit = signatures / norms
    return unit @ unit.T
