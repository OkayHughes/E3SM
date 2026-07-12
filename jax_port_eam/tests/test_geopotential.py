"""Tier-1 golden replay + Tier-0 properties for eam_jax.geopotential."""

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
GOLDEN = ROOT / "golden" / "geopotential_golden.npz"

from eam_jax.geopotential import geopotential_t  # noqa: E402


@pytest.fixture(scope="module")
def gold():
    if not GOLDEN.exists():
        pytest.skip("golden archive missing (run harness in container)")
    return np.load(GOLDEN)


def test_golden_replay(gold):
    zi, zm = geopotential_t(gold["piln"], gold["pmln"], gold["pint"],
                            gold["pmid"], gold["pdel"], gold["rpdel"],
                            gold["t"], gold["q"], 287.042, 9.80616,
                            0.60779262)
    np.testing.assert_allclose(np.asarray(zi), gold["zi_t"], rtol=1e-13,
                               atol=1e-10)
    np.testing.assert_allclose(np.asarray(zm), gold["zm_t"], rtol=1e-13)


def test_properties(gold):
    zi, zm = geopotential_t(gold["piln"], gold["pmln"], gold["pint"],
                            gold["pmid"], gold["pdel"], gold["rpdel"],
                            gold["t"], gold["q"], 287.042, 9.80616,
                            0.60779262)
    zi, zm = np.asarray(zi), np.asarray(zm)
    # surface interface is exactly zero; heights decrease downward
    assert (zi[:, -1] == 0.0).all()
    assert (np.diff(zi, axis=1) < 0).all()
    # midpoints sit between their bounding interfaces
    assert (zm > zi[:, 1:]).all() and (zm < zi[:, :-1]).all()
    # scale height sanity: column top of a 72L ~1000-225.5Pa grid
    assert 15e3 < zi[:, 0].min() and zi[:, 0].max() < 60e3
