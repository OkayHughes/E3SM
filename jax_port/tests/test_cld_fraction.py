"""Tier-0 tests for scream_jax.cld_fraction.

Cross-checks the JAX kernel against EAMxx's OWN NumPy reference
implementation (components/eamxx/src/physics/cld_fraction/
cld_fraction_numpy.py) — upstream code, imported directly, so this is a
genuine (if easy) golden comparison. The in-situ Tier-2 comparison against
the C++ runs inside EAMxx (see TEST_HARNESS_DESIGN.md).
"""

import importlib.util
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "jax_port"))

from scream_jax.cld_fraction import cld_fraction_main  # noqa: E402
from scream_jax.adapters.eamxx import cld_fraction_jax  # noqa: E402


def _load_upstream_numpy_ref():
    path = (REPO / "components" / "eamxx" / "src" / "physics"
            / "cld_fraction" / "cld_fraction_numpy.py")
    spec = importlib.util.spec_from_file_location("cld_fraction_numpy", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _random_inputs(seed, shape=(218, 72)):
    rng = np.random.default_rng(seed)
    qi = 10.0 ** rng.uniform(-9.0, -3.0, shape)          # spans both thresholds
    qi[rng.uniform(size=shape) < 0.3] = 0.0              # plenty of clear sky
    liq = np.clip(rng.uniform(-0.3, 1.0, shape), 0.0, 1.0)
    return qi, liq


def test_matches_upstream_numpy_reference():
    ref = _load_upstream_numpy_ref()
    qi, liq = _random_inputs(0)
    thr, thr4out = 1e-12, 1e-5  # EAMxx defaults (ice_cloud_threshold / _for_analysis)

    ice_r = np.empty_like(qi); tot_r = np.empty_like(qi)
    ice4_r = np.empty_like(qi); tot4_r = np.empty_like(qi)
    ref.main(thr, thr4out, qi, liq, ice_r, tot_r, ice4_r, tot4_r)

    ice_j, tot_j, ice4_j, tot4_j = (np.asarray(a) for a in
                                    cld_fraction_main(thr, thr4out, qi, liq))

    np.testing.assert_array_equal(ice_j, ice_r)
    np.testing.assert_array_equal(tot_j, tot_r)
    np.testing.assert_array_equal(ice4_j, ice4_r)
    np.testing.assert_array_equal(tot4_j, tot4_r)


def test_adapter_in_place_semantics():
    # The adapter must fill preallocated arrays exactly like the reference.
    ref = _load_upstream_numpy_ref()
    qi, liq = _random_inputs(1, shape=(8, 12))
    thr, thr4out = 1e-12, 1e-5

    outs_a = [np.full_like(qi, -99.0) for _ in range(4)]
    outs_r = [np.full_like(qi, -99.0) for _ in range(4)]
    cld_fraction_jax.main(thr, thr4out, qi, liq, *outs_a)
    ref.main(thr, thr4out, qi, liq, *outs_r)

    for a, r in zip(outs_a, outs_r):
        np.testing.assert_array_equal(a, r)
        assert not np.any(a == -99.0)


def test_threshold_is_strict_inequality():
    # qi exactly at the threshold is NOT cloudy (C++ uses qi > threshold).
    thr = 1e-12
    qi = np.array([[thr, np.nextafter(thr, 1.0), 0.0]])
    liq = np.zeros_like(qi)
    ice, tot, _, _ = (np.asarray(a) for a in cld_fraction_main(thr, 1e-5, qi, liq))
    np.testing.assert_array_equal(ice, [[0.0, 1.0, 0.0]])
    np.testing.assert_array_equal(tot, ice)
