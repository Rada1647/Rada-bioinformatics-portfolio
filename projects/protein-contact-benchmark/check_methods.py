"""Independent small mathematical/implementation checks, without protein data.

Run: python check_methods.py
Checks use dense linear algebra and explicit endpoint-triple enumeration,
independently of the production banded resolvent and cubic formula.
"""

from __future__ import annotations

import itertools
import json

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
from scipy.ndimage import gaussian_filter1d, median_filter, uniform_filter1d

from methods import candidates, predict, predict_all_for_validation, _power


def path_matrix(n):
    A = np.eye(n)
    for i in range(n - 1):
        A[i, i] -= 0.25
        A[i + 1, i + 1] -= 0.25
        A[i, i + 1] = A[i + 1, i] = 0.25
    return A


def rejected(call):
    try:
        call()
    except ValueError:
        return
    raise AssertionError("invalid input was accepted")


def run_checks():
    checks = []
    rng = np.random.default_rng(3322026)
    for n in (1, 2, 3, 9):
        Z = rng.integers(0, 2, (n, 5)).astype(np.float64)
        A = path_matrix(n)
        assert_allclose(A.sum(axis=1), 1, atol=1e-15)
        assert_array_equal(A, A.T)
        for t in (1, 8, 64):
            G = np.linalg.inv(np.eye(n) + t * (np.eye(n) - A))
            H = predict(Z, "linear", {"t": t})
            assert_allclose(H, G @ Z, atol=3e-14, rtol=3e-14)
            assert_allclose(H.mean(axis=0), Z.mean(axis=0), atol=3e-14)
    checks.append("banded resolvent agrees with independently constructed dense inverse; uniform occupancy preserved")

    Z = np.array([[0., 1.], [1., 0.], [1., 1.]])
    n = len(Z)
    for t in (1, 8, 64):
        A = path_matrix(n)
        G = np.linalg.inv(np.eye(n) + t * (np.eye(n) - A))
        expected = np.zeros_like(Z)
        for i in range(n):
            for a, b, c in itertools.product(range(n), repeat=3):
                expected[i] += G[i, a] * G[i, b] * G[i, c] * np.median(Z[[a, b, c]], axis=0)
        assert_allclose(predict(Z, "paper332", {"t": t}), expected, atol=4e-14, rtol=4e-14)
    checks.append("paper construction agrees with exhaustive expected median of three independent endpoints")

    grid = candidates()
    assert sum(map(len, grid.values())) == 120
    assert len(grid["power_full"]) == 63 and len(grid["power_matched"]) == 7
    assert {"t": 8, "alpha": 3.0} in grid["power_matched"]
    mutable = candidates()
    mutable["linear"][0]["t"] = 999
    assert candidates()["linear"][0]["t"] == 1
    Z = rng.integers(0, 2, (17, 4)).astype(np.float64)
    original = Z.copy()
    seen = set()
    key = lambda m, p: (m, json.dumps(p, sort_keys=True))
    for method, params, P in predict_all_for_validation(Z, noise_level=0.1):
        assert P.dtype == np.float64 and P.shape == Z.shape
        assert np.all(np.isfinite(P)) and P.min() >= 0 and P.max() <= 1
        assert_allclose(P, predict(Z, method, params, noise_level=0.1), atol=0, rtol=0)
        assert key(method, params) not in seen
        seen.add(key(method, params))
    assert seen == {key(m, p) for m, pp in grid.items() for p in pp}
    assert_array_equal(Z, original)
    assert sum(1 for _ in predict_all_for_validation(Z)) == 113
    checks.append("all 120 candidate outputs match scalar API exactly; 113 candidates without oracle; input unchanged")

    for method, params_list in grid.items():
        for params in params_list:
            for n in (1, 11):
                C = np.column_stack((np.zeros(n), np.ones(n)))
                P = predict(C, method, params, noise_level=0.1)
                assert_allclose(P, C, atol=3e-14, rtol=0)
    checks.append("constants and singleton trajectories preserved for every candidate")

    # A fixed left segment must not depend on values from an independent run.
    left = np.array([[0.], [0.], [1.], [0.]])
    for method, params_list in grid.items():
        p = params_list[-1]
        baseline = predict(left, method, p, noise_level=0.1)
        for right in (np.zeros((7, 1)), np.ones((7, 1))):
            outputs = [predict(segment, method, p, noise_level=0.1) for segment in (left, right)]
            assert_array_equal(outputs[0], baseline)
    checks.append("per-trajectory calls have no shared state or cross-sequence contamination")

    Z = np.array([[1., 0.], [0., 1.], [0., 0.], [0., 1.], [0., 0.]])
    assert_allclose(predict(Z, "moving_average", {"half_width": 2}),
                    uniform_filter1d(Z, 5, axis=0, mode="reflect"), atol=0)
    assert_allclose(predict(Z, "gaussian", {"sigma": 1.0}),
                    gaussian_filter1d(Z, 1.0, axis=0, mode="reflect", truncate=4.0), atol=0)
    assert_array_equal(predict(Z, "median", {"width": 3}),
                       median_filter(Z, size=(3, 1), mode="reflect"))
    # Independent explicit reflected padding checks the boundary convention.
    padding = np.pad(Z, ((2, 2), (0, 0)), mode="symmetric")
    explicit = np.stack([padding[i:i + 5].mean(axis=0) for i in range(len(Z))])
    assert_allclose(predict(Z, "moving_average", {"half_width": 2}), explicit, atol=2e-16)
    padding = np.pad(Z, ((1, 1), (0, 0)), mode="symmetric")
    explicit = np.stack([np.median(padding[i:i + 3], axis=0) for i in range(len(Z))])
    assert_array_equal(predict(Z, "median", {"width": 3}), explicit)
    checks.append("filters match specified reflect boundaries and do not mix features")

    H = np.array([0.0, 1e-12, .01, .25, .499999, .5, .500001, .75, .99, 1 - 1e-12, 1.0])
    cubic = H * H * (3 - 2 * H)
    assert_array_equal(cubic >= .5, H >= .5)
    assert np.all(np.diff(cubic) >= 0)
    for alpha in (.5, .75, 1, 1.25, 1.5, 2, 3, 4, 8):
        P = _power(H, alpha)
        direct = H ** alpha / (H ** alpha + (1 - H) ** alpha)
        assert_allclose(P, direct, atol=3e-16, rtol=2e-14)
        assert P[0] == 0 and P[-1] == 1 and P[5] == .5
        assert_array_equal(P >= .5, H >= .5)
        assert np.all(np.diff(P) >= 0)
    checks.append("cubic and positive-power transforms preserve ranking, threshold, and exact ties")

    Z = np.array([[0.], [1.]])
    H = predict(Z, "linear", {"t": 4})
    assert_allclose(predict(Z, "oracle_debiased_linear", {"t": 4}, noise_level=.1),
                    np.clip((H - .1) / .8, 0, 1), atol=0)
    assert_array_equal(predict(Z, "hard_linear", {"t": 4}), (H >= .5).astype(float))
    checks.append("known-noise correction and hard threshold agree with declared formulas")

    for bad in ([], [0, 1], [[.5]], [[np.nan]], [[np.inf]], [[1j]], np.empty((2, 0))):
        rejected(lambda bad=bad: predict(bad, "paper332", {"t": 1}))
    for bad in (0, -1, 1.1, True, float("nan")):
        rejected(lambda bad=bad: predict(Z, "linear", {"t": bad}))
    for bad in (None, -.1, .5, 1, float("nan"), True):
        rejected(lambda bad=bad: predict(Z, "oracle_debiased_linear", {"t": 1}, noise_level=bad))
    rejected(lambda: predict(Z, "median", {"width": 2}))
    rejected(lambda: predict(Z, "power_full", {"t": 1, "alpha": 0}))
    rejected(lambda: predict(Z, "linear", {"t": 1, "ignored": 7}))
    rejected(lambda: predict(Z, "unknown", {}))
    checks.append("malformed, continuous, nonfinite inputs and invalid parameters rejected")
    return {"status": "passed", "checks": checks, "candidate_count_with_oracle": 120,
            "candidate_count_without_oracle": 113, "real_protein_data_used": False}


if __name__ == "__main__":
    print(json.dumps(run_checks(), indent=2))
