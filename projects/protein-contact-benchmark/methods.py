"""Predeclared binary-contact smoothing methods for the #332 extension.

The lazy reflecting path has interior diagonal 1/2 and adjacent weights 1/4;
endpoints have diagonal 3/4. A singleton path is the identity. For binary
coordinates the manuscript construction is exactly, in real arithmetic,
H = [I + t(I-A)]**(-1) Z, Y = H**2 * (3 - 2*H). This is the expectation of
the coordinatewise median of three iid geometric-walk endpoints, not an
arbitrary transformation of a continuous signal.

The resolvent implementation is adapted from ../adk_contact332/smoothing.py.
Mathematical source, pinned to the release used in the original benchmark:
https://github.com/openai/math/blob/adc7f1241b42e322a6451854ab7e4b4c146bf78a/
preprints/Metric-Markov-Cotype-Two-of-l1-October-5-2026/l1-markov-cotype.pdf

Each call handles ONE uninterrupted trajectory or split; callers must never
concatenate independent trajectories before smoothing. Axis 0 is time.
All methods are retrospective: future and past observations can contribute.
All numerical calculations and outputs use float64. The statistical Brier
benchmark is not a test of the paper's squared-whole-vector-L1 inequality.
"""

from __future__ import annotations

from numbers import Integral, Real

import numpy as np
from scipy.linalg import solve_banded
from scipy.ndimage import gaussian_filter1d, median_filter, uniform_filter1d
from scipy.special import expit, logit


T_GRID = (1, 2, 4, 8, 16, 32, 64)
ALPHA_GRID = (0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0, 4.0, 8.0)
# Alpha=3 was selected on the prior ADK development experiment; only the
# smoothing parameter is tuned here, giving the same seven-candidate budget
# as the paper construction. This grid is fixed before new-protein testing.
MATCHED_POWER_GRID = tuple((t, 3.0) for t in T_GRID)
SIGMA_GRID = (0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0)
MEDIAN_WIDTH_GRID = (3, 5, 9, 17, 33, 65, 129)
_TOL = 1e-11


def candidates():
    """Return fresh JSON-serializable parameter grids in deterministic order.

    There are 120 configurations when the oracle is eligible, or 113 without
    it. Candidate order resolves exact validation-score ties within families.
    Power-full order is ascending t, then ascending alpha.
    """
    return {
        "raw": [{}],
        "linear": [{"t": t} for t in T_GRID],
        "paper332": [{"t": t} for t in T_GRID],
        "oracle_debiased_linear": [{"t": t} for t in T_GRID],
        "hard_linear": [{"t": t} for t in T_GRID],
        "moving_average": [{"half_width": w} for w in T_GRID],
        "gaussian": [{"sigma": s} for s in SIGMA_GRID],
        "median": [{"width": w} for w in MEDIAN_WIDTH_GRID],
        "power_full": [{"t": t, "alpha": a} for t in T_GRID for a in ALPHA_GRID],
        "power_matched": [{"t": t, "alpha": a} for t, a in MATCHED_POWER_GRID],
    }


def _binary_data(Z):
    if np.iscomplexobj(Z):
        raise ValueError("Z must be real binary data")
    Z = np.asarray(Z, dtype=np.float64)
    if Z.ndim != 2 or min(Z.shape) < 1 or not np.all(np.isfinite(Z)):
        raise ValueError("Z must be a nonempty finite two-dimensional array")
    if not np.all((Z == 0.0) | (Z == 1.0)):
        raise ValueError("Z must contain only binary values 0 and 1")
    return Z


def _positive_integer(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral) or value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def _positive_real(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real) or not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be a finite positive number")
    return float(value)


def _epsilon(noise_level):
    if isinstance(noise_level, (bool, np.bool_)) or not isinstance(noise_level, Real) or not np.isfinite(noise_level) or not 0 <= noise_level < 0.5:
        raise ValueError("oracle noise_level must be finite and satisfy 0 <= epsilon < 0.5")
    return float(noise_level)


def _resolvent(Z, t):
    """Tridiagonal solve for an already validated binary float64 input."""
    n = Z.shape[0]
    if n == 1:
        return Z.copy()
    ab = np.zeros((3, n), dtype=np.float64)
    ab[0, 1:] = -0.25 * t
    ab[2, :-1] = -0.25 * t
    ab[1] = 1.0 + 0.5 * t
    ab[1, [0, -1]] = 1.0 + 0.25 * t
    H = solve_banded((1, 1), ab, Z, check_finite=False)
    if np.min(H) < -_TOL or np.max(H) > 1 + _TOL:
        raise ArithmeticError("resolvent probabilities exceeded numerical tolerance")
    H = np.clip(H, 0.0, 1.0)
    # The resolvent fixes constant columns exactly. Restore their known
    # endpoints so alpha<1 sharpening does not amplify solver roundoff near
    # 1 into an artificial ~sqrt(machine-epsilon) change in a constant input.
    occupied = np.sum(Z, axis=0)
    H[:, occupied == 0] = 0.0
    H[:, occupied == n] = 1.0
    return H


def _cubic(H):
    return H * H * (3.0 - 2.0 * H)


def _power(H, alpha):
    """Stable h**alpha / (h**alpha + (1-h)**alpha), including endpoints."""
    if alpha == 1.0:
        return H.copy()
    return expit(alpha * logit(H))


def _parameters(method, params):
    """Validate the form, allowing off-grid values for mathematical checks."""
    keys = {
        "raw": (), "linear": ("t",), "paper332": ("t",),
        "oracle_debiased_linear": ("t",), "hard_linear": ("t",),
        "moving_average": ("half_width",), "gaussian": ("sigma",),
        "median": ("width",), "power_full": ("t", "alpha"),
        "power_matched": ("t", "alpha"),
    }
    if method not in keys:
        raise ValueError(f"unknown method: {method}")
    if not isinstance(params, dict) or set(params) != set(keys[method]):
        raise ValueError(f"{method} requires exactly these parameter keys: {keys[method]}")
    out = {}
    for name, value in params.items():
        out[name] = (_positive_integer(value, name) if name in ("t", "half_width", "width")
                     else _positive_real(value, name))
    if "width" in out and out["width"] % 2 != 1:
        raise ValueError("median width must be odd")
    return out


def predict(Z, method, params, noise_level=None):
    """Return same-shape float64 probabilities for one binary trajectory.

    Moving average, Gaussian, and median filters use SciPy's ``reflect``
    (half-sample symmetric) extension at each trajectory boundary. Gaussian
    kernels use truncate=4.0, order=0. Median filtering uses size=(width, 1),
    so features never mix. Hard calls use >= 0.5, including exact ties.
    The oracle baseline requires the known symmetric bit-flip probability;
    it is not applicable to coordinate-jitter corruption. No metadata allows
    this function to distinguish corruption mechanisms: callers enforce that.
    """
    Z = _binary_data(Z)
    params = _parameters(method, params)
    if method == "raw":
        return Z.copy()
    if method == "moving_average":
        return np.clip(uniform_filter1d(Z, size=2 * params["half_width"] + 1,
                                       axis=0, mode="reflect"), 0.0, 1.0)
    if method == "gaussian":
        return np.clip(gaussian_filter1d(Z, sigma=params["sigma"], axis=0,
                                        order=0, mode="reflect", truncate=4.0), 0.0, 1.0)
    if method == "median":
        return median_filter(Z, size=(params["width"], 1), mode="reflect")
    if method == "oracle_debiased_linear":
        epsilon = _epsilon(noise_level)
    H = _resolvent(Z, params["t"])
    if method == "linear":
        return H
    if method == "paper332":
        return _cubic(H)
    if method == "oracle_debiased_linear":
        return np.clip((H - epsilon) / (1.0 - 2.0 * epsilon), 0.0, 1.0)
    if method == "hard_linear":
        return (H >= 0.5).astype(np.float64)
    return _power(H, params["alpha"])


def predict_all_for_validation(Z, noise_level=None):
    """Yield ``(method, params, prediction)`` for every eligible candidate.

    Outputs are grouped by t to share each of seven resolvent solves. This
    generator order differs from candidates(); consumers must use method and
    params as identities and break validation ties with candidates() order.
    The oracle family is omitted when noise_level is None. At most one H is
    retained by the generator, rather than caching seven full trajectories.
    Consumers must not mutate yielded predictions before resuming iteration.
    """
    Z = _binary_data(Z)
    epsilon = None if noise_level is None else _epsilon(noise_level)
    yield "raw", {}, Z.copy()
    matched = dict(MATCHED_POWER_GRID)
    for t in T_GRID:
        H = _resolvent(Z, t)
        yield "linear", {"t": t}, H
        yield "paper332", {"t": t}, _cubic(H)
        if epsilon is not None:
            yield "oracle_debiased_linear", {"t": t}, np.clip((H - epsilon) / (1 - 2 * epsilon), 0.0, 1.0)
        yield "hard_linear", {"t": t}, (H >= 0.5).astype(np.float64)
        for alpha in ALPHA_GRID:
            yield "power_full", {"t": t, "alpha": alpha}, _power(H, alpha)
        yield "power_matched", {"t": t, "alpha": matched[t]}, _power(H, matched[t])
    for method in ("moving_average", "gaussian", "median"):
        for params in candidates()[method]:
            yield method, params, predict(Z, method, params)
