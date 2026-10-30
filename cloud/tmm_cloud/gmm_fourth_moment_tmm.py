"""Freeze a fitted GMM and identify one t degree of freedom per component.

01_IDEA/derivation.md section 7, equations (33)-(35). Input is the same joint
state coordinates used by the saved GMM. No t-EM, repartitioning or covariance
refitting is performed. This module is an addition to the frozen cloud package.
"""

import copy
import math
import time


def df_from_radial_fourth(beta, dimension):
    """Covariance-standardized Mardia fourth moment -> df and excess ratio.

The Gaussian boundary is explicit. A below-Gaussian fourth moment has no
finite-df solution in this family; infinity is a boundary model, not a match.
"""
    if isinstance(dimension, bool) or not isinstance(dimension, int) or dimension < 1:
        raise ValueError("dimension must be a positive integer")
    if not math.isfinite(beta) or beta < 0:
        raise ValueError("radial fourth moment must be finite and nonnegative")
    gaussian = dimension * (dimension + 2)
    excess = beta / gaussian - 1.0
    if excess <= 0:
        return math.inf, excess
    nu = 4.0 + 2.0 / excess
    if not math.isfinite(nu) or not nu > 4.0:
        raise ArithmeticError("positive fourth-moment excess cannot be represented as finite df > 4")
    return nu, excess


def covariance_to_scatter_factor(nu):
    if nu == math.inf:
        return 1.0
    if not math.isfinite(nu) or not nu > 4:
        raise ValueError("fourth-moment t calibration requires finite df > 4 or positive infinity")
    return 1.0 - 2.0 / nu


def fit_gmm_fourth_moment_tmm(Z, gmm, chunk_size=1024):
    """Use actual modeling states and frozen Gaussian posterior probabilities.

Requires the existing cloud/tmm_cloud/mixture.py from source commit 32bc9ae.
The result is the same Mixture class and can use the original density and
propagation functions. Gaussian components have df=inf, as in the source API.
"""
    import numpy as np
    from . import mixture as mx

    started = time.perf_counter()
    Z = np.asarray(Z, dtype=float)
    if Z.ndim != 2 or not np.isfinite(Z).all() or len(Z) == 0:
        raise ValueError("modeling coordinates must be a finite nonempty matrix")
    n, d = Z.shape
    if not np.isinf(gmm.df).all() or not (gmm.df > 0).all():
        raise ValueError("source model must be a Gaussian mixture")
    if gmm.means.shape != (gmm.K, d) or gmm.scatter.shape != (gmm.K, d, d):
        raise ValueError("GMM dimensions do not match the original coordinates")
    if not (gmm.weights > 0).all() or not math.isclose(float(gmm.weights.sum()), 1.0, abs_tol=1e-10):
        raise ValueError("source GMM weights must be positive and sum to one")
    if not isinstance(chunk_size, int) or chunk_size < 1:
        raise ValueError("chunk_size must be a positive integer")

    counts = np.zeros(gmm.K)
    squared_membership = np.zeros(gmm.K)
    fourth_sum = np.zeros(gmm.K)
    for start in range(0, n, chunk_size):
        # Source Gaussian Mahalanobis distance is already a squared distance.
        logs, distance = mx.mixture_terms(Z[start:start + chunk_size], gmm)
        membership = np.exp(logs - mx.logsumexp_rows(logs)[:, None])
        counts += membership.sum(axis=0)
        squared_membership += (membership * membership).sum(axis=0)
        fourth_sum += (membership * distance * distance).sum(axis=0)
    if not (counts > 0).all() or not np.isfinite(fourth_sum).all():
        raise ValueError("at least one component lacks finite fourth-moment evidence")

    model = copy.deepcopy(gmm)
    model.family = "tmm"
    model.df = np.empty(gmm.K)
    model.scatter = gmm.scatter.copy()
    components = []
    gaussian_beta = d * (d + 2)
    for k in range(gmm.K):
        beta = float(fourth_sum[k] / counts[k])
        nu, excess = df_from_radial_fourth(beta, d)
        scale = covariance_to_scatter_factor(nu)
        model.df[k] = nu
        model.scatter[k] *= scale
        reconstructed = model.scatter[k] if math.isinf(nu) else model.scatter[k] * (nu / (nu - 2))
        if not np.allclose(reconstructed, gmm.scatter[k], rtol=1e-12, atol=1e-14):
            raise ArithmeticError("component covariance changed during the t conversion")
        theoretical_beta = gaussian_beta if math.isinf(nu) else gaussian_beta * (nu - 2) / (nu - 4)
        components.append(dict(component=k, dimension=d, posterior_count=float(counts[k]),
                               kish_effective_count=float(counts[k] ** 2 / squared_membership[k]),
                               radial_fourth=beta, gaussian_radial_fourth=gaussian_beta,
                               radial_excess=excess, df=None if math.isinf(nu) else nu,
                               df_kind="Gaussian limit" if math.isinf(nu) else "finite t",
                               scatter_factor=scale, theoretical_radial_fourth=theoretical_beta,
                               fourth_moment_residual=theoretical_beta - beta,
                               fourth_moment_matched=excess >= 0))
    if not np.array_equal(model.weights, gmm.weights) or not np.array_equal(model.means, gmm.means):
        raise ArithmeticError("source weights or means changed")
    model.audit = dict(method="GMM-fixed fourth-moment TMM", K=gmm.K, dimension=d,
                       status="closed-form", accepted=True, modeling_scenarios=n,
                       posterior_source="frozen source GMM; no t posterior update",
                       fourth_statistic="posterior-weighted covariance-standardized radial fourth moment",
                       df_rule="4+2/excess if excess>0; Gaussian limit otherwise",
                       covariance_preserved=True, weights_preserved=True, means_preserved=True,
                       components=components, seconds=time.perf_counter() - started,
                       evidence="calibration only; AC distribution accuracy requires independent comparison")
    return model
