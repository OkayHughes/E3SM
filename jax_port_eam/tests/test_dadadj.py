"""Tier-1 golden replay + Tier-0 properties for eam_jax.dadadj."""

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
GOLDEN = ROOT / "golden" / "dadadj_golden.npz"

from eam_jax.dadadj import dadadj  # noqa: E402


@pytest.fixture(scope="module")
def gold():
    if not GOLDEN.exists():
        pytest.skip("golden archive missing (run harness in container)")
    return np.load(GOLDEN)


@pytest.mark.parametrize("nlvdry", [3, 8])
def test_golden_replay(gold, nlvdry):
    t, q, conv = dadadj(gold["pmid"], gold["pint"], gold["pdel"],
                        gold[f"t_in_{nlvdry}"], gold[f"q_in_{nlvdry}"],
                        nlvdry=nlvdry)
    assert bool(np.asarray(conv).all())
    np.testing.assert_allclose(np.asarray(t), gold[f"t_out_{nlvdry}"],
                               rtol=1e-13)
    np.testing.assert_allclose(np.asarray(q), gold[f"q_out_{nlvdry}"],
                               rtol=1e-13, atol=1e-300)


def test_stable_profile_untouched(gold):
    # stable columns (0..3 in the generator) must be bit-identical
    t, q, conv = dadadj(gold["pmid"], gold["pint"], gold["pdel"],
                        gold["t_in_3"], gold["q_in_3"], nlvdry=3)
    np.testing.assert_array_equal(np.asarray(t)[:4], gold["t_in_3"][:4])
    np.testing.assert_array_equal(np.asarray(q)[:4], gold["q_in_3"][:4])


def test_conservation_and_stability(gold):
    """Enthalpy-like sum(c_p-weighted...) — dadadj conserves
    pdel-weighted dry static energy proxy sum(t*pdel*(c2 relation)) is
    intricate; the robust invariants are: q mass conserved, adjusted
    region no longer superadiabatic beyond tolerance."""
    from eam_jax.constants import CAPPA
    nlvdry = 8
    pmid, pint, pdel = gold["pmid"], gold["pint"], gold["pdel"]
    t0, q0 = gold[f"t_in_{nlvdry}"], gold[f"q_in_{nlvdry}"]
    t, q, conv = dadadj(pmid, pint, pdel, t0, q0, nlvdry=nlvdry)
    t, q = np.asarray(t), np.asarray(q)
    # water mass conservation over the adjusted band (pairs mix q)
    band = slice(0, nlvdry + 1)
    np.testing.assert_allclose((q[:, band] * pdel[:, band]).sum(1),
                               (q0[:, band] * pdel[:, band]).sum(1),
                               rtol=1e-12)
    # post-adjustment: no pair superadiabatic beyond the final zeps
    k = np.arange(nlvdry)
    gammad = CAPPA * 0.5 * (t[:, k + 1] + t[:, k]) / pint[:, k + 1]
    dtdp = (t[:, k + 1] - t[:, k]) / (pmid[:, k + 1] - pmid[:, k])
    assert (dtdp <= gammad + 2e-4).all()
