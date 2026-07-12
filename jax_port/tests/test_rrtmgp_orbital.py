"""Orbital mechanics (shr_orb) vs Fortran reference values, and trcmix
sanity."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scream_jax.rrtmgp import orbital as orb  # noqa: E402


def test_shr_orb_vs_fortran():
    # Reference values from share/util/shr_orb_mod.F90 compiled in the
    # dev container (year 2021, calday 285.520833)
    eccen, obliq, mvelp, obliqr, lambm0, mvelpp = orb.shr_orb_params(2021)
    assert abs(eccen - 1.66951276950526098e-02) < 1e-16
    assert abs(mvelpp - 4.94373707898766490e+00) < 1e-13
    assert abs(lambm0 - -3.24061841596024122e-02) < 1e-15
    assert abs(obliqr - 4.09053448128890051e-01) < 1e-14

    delta, eccf = orb.shr_orb_decl(285.520833, eccen, mvelpp, lambm0, obliqr)
    assert abs(delta - -1.26227807149654370e-01) < 1e-14
    assert abs(eccf - 1.00358634673266378e+00) < 1e-14

    refs = [(0.1, 0.3, 8.84290934066529810e-01, 8.54785866181254517e-01),
            (-0.8, 2.9, -5.96623358604308152e-01, 0.0),
            (1.2, -1.0, 1.14708336815450646e-01, 1.31989765436686640e-01)]
    for la, lo, r0, r1 in refs:
        c0 = float(orb.shr_orb_cosz(285.520833, la, lo, delta, 0.0))
        c1 = float(orb.shr_orb_cosz(285.520833, la, lo, delta, 1800.0))
        assert abs(c0 - r0) < 1e-14
        assert abs(c1 - r1) < 1e-14


def test_trcmix_profiles():
    lat = np.array([0.0, 45.0, -80.0])
    pmid = np.linspace(100.0, 1000e2, 40)[None, :] * np.ones((3, 1))
    for gas in ("co2", "o2", "ch4", "n2o", "cfc11", "cfc12"):
        q = orb.trcmix(gas, lat, pmid, 400e-6, 320e-9, 1.8e-6, 7.7e-10,
                       5.3e-10)
        assert np.isfinite(q).all() and (q > 0).all()
        # tropospheric value constant, decreasing into the stratosphere
        if gas not in ("co2", "o2"):
            assert q[0, -1] >= q[0, 0]
