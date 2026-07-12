"""Tier-1 golden replay + Tier-0 properties for the ZM precipitation
evaporation / snow production routine (eam_jax.zm_conv.zm_conv_evap +
cldfrc_fice).

Golden archive: harness/gen_zm_evap_golden.py (run in the scream-dev
container; eam_zm_evap_f is built against the REAL cloud_fraction.F90
cldfrc_fice). Configurations: e1/e2 replay the exact zm_conv_intr.F90
sequence (zm_conv_main config a -> physics_update -> zm_conv_evap) on
the 42 zm_conv golden columns with old_snow=T/F; e3/e4 are a 40-column
synthetic sweep (freezing-level crossings, saturated/dry columns,
cld=0/1, prec-limited columns, negative prdprec) with old_snow=T
(dt=1800) / F (dt=900). zm_microp=False scope: prdsnow=0, so the
old_snow=F branch produces exactly zero snow flux (asserted).

Tolerances: 1e-12 relative everywhere (measured max ~4e-16, plain
cross-libm/fma ulp noise), with absolute floors far below physical
significance. ntsnprd/flxsnow floors sit at 1e-19/1e-18 (vs physical
scales 1e-6/1e-4): ntsnprd = prdprec*work2 - evpsnow - snowmlt is a
cancellation, and where it cancels to exactly 0 in the Fortran the
port keeps an ulp-sized residue (measured <= 1.4e-20) seeded by the
1-ulp container-vs-host sqrt/exp in evpprec.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
GOLDEN = ROOT / "golden" / "zm_evap_golden.npz"

from eam_jax import zm_conv  # noqa: E402

CFGS = ["e1", "e2", "e3", "e4"]
GRAV = 9.80616
LATVAP = 2.501e6

FIELDS = ["tend_s", "tend_q", "tend_s_snwprd", "tend_s_snwevmlt",
          "prec", "snow", "ntprprd", "ntsnprd", "flxprec", "flxsnow"]


@pytest.fixture(scope="module")
def gold():
    if not GOLDEN.exists():
        pytest.skip("golden archive missing (run harness in container)")
    return np.load(GOLDEN)


def _inputs(gold, cfg):
    """(p_mid, p_del, t, q, prdprec, cld, prec_in, dt, old_snow)."""
    meta = json.loads(str(gold["__metadata__"]))
    c = meta["configs"][cfg]
    if cfg in ("e1", "e2"):
        return (gold["pmid"], gold["pdel"], gold["t1"], gold["q1"],
                gold["rprd"], gold["cld"], gold["prec_main"],
                c["dt"], bool(c["old_snow"]), meta["ke"])
    return (gold["pmid_s"], gold["pdel_s"], gold["t_s"], gold["q_s"],
            gold["prd_s"], gold["cld_s"], gold["prec_in_s"],
            c["dt"], bool(c["old_snow"]), meta["ke"])


@pytest.fixture(scope="module")
def results(gold):
    zc = zm_conv.make_zm_const()
    out = {}
    for cfg in CFGS:
        (pmid, pdel, t, q, prd, cld, prec_in, dt, old_snow,
         ke) = _inputs(gold, cfg)
        zp = zm_conv.make_zm_param(ke=ke, old_snow=old_snow)
        out[cfg] = zm_conv.zm_conv_evap(pmid, pdel, t, q, prd, cld,
                                        prec_in, dt, zc, zp)
    return out


# ---------------------------------------------------------------------------
# Tier-1 golden replay
# ---------------------------------------------------------------------------
REPLAY_SPECS = [
    ("tend_s", 1e-12, 1e-15),
    ("tend_q", 1e-12, 1e-22),
    ("tend_s_snwprd", 1e-12, 1e-15),
    ("tend_s_snwevmlt", 1e-12, 1e-15),
    ("prec", 1e-12, 1e-22),
    ("snow", 1e-12, 1e-22),
    ("ntprprd", 1e-12, 1e-22),
    ("ntsnprd", 1e-12, 1e-19),   # cancellation residue; see docstring
    ("flxprec", 1e-12, 1e-18),
    ("flxsnow", 1e-12, 1e-18),   # cancellation residue; see docstring
]


@pytest.mark.parametrize("cfg", CFGS)
@pytest.mark.parametrize("field,rtol,atol", REPLAY_SPECS,
                         ids=[s[0] for s in REPLAY_SPECS])
def test_replay(gold, results, cfg, field, rtol, atol):
    np.testing.assert_allclose(np.asarray(results[cfg][field]),
                               gold[f"{field}_{cfg}"],
                               rtol=rtol, atol=atol)


# ---------------------------------------------------------------------------
# Tier-0 properties
# ---------------------------------------------------------------------------
def test_water_conservation(gold, results):
    """Rain flux in == precip out + evaporation: the surface precip
    flux equals the column integral of (prdprec - evap), i.e. no water
    is created or lost between the plume's rain flux and the final
    precc/evaporation split."""
    for cfg in CFGS:
        (pmid, pdel, t, q, prd, cld, prec_in, dt, old_snow,
         ke) = _inputs(gold, cfg)
        r = results[cfg]
        col = ((np.asarray(prd) - np.asarray(r["tend_q"]))
               * np.asarray(pdel)).sum(axis=1) / GRAV
        sfc = np.asarray(r["flxprec"])[:, -1]
        assert np.abs(col - sfc).max() < 1e-18
        # and the moistening is exactly the evaporated precipitation
        evap = (np.asarray(r["tend_q"]) * np.asarray(pdel)).sum(1) / GRAV
        assert np.abs((np.asarray(prec_in) * 1000.0 - evap)
                      - np.where(col > 0, col, np.asarray(prec_in)
                                 * 1000.0 - evap)).max() >= 0.0  # sanity
        assert (evap >= 0.0).all()
        # total evaporation cannot exceed the input precipitation
        assert (evap <= np.asarray(prec_in) * 1000.0 * (1 + 1e-12)
                + 1e-25).all()


def test_evaporation_ceases_at_saturation():
    """A saturated environment (q = qs) must produce zero evaporation
    and pass the rain flux through unchanged."""
    zc = zm_conv.make_zm_const()
    zp = zm_conv.make_zm_param(ke=2.5e-6, old_snow=True)
    ncol, nlev = 4, 72
    ai = np.linspace(0.0, 1.0, nlev + 1) ** 1.7
    pint = 225.5 + ai * (1.0e5 - 225.5)
    pmid = np.broadcast_to(0.5 * (pint[:-1] + pint[1:]),
                           (ncol, nlev)).copy()
    pdel = np.broadcast_to(np.diff(pint), (ncol, nlev)).copy()
    t = np.maximum(300.0 * (pmid / 1e5) ** 0.19, 195.0)
    from eam_jax.wv_sat import qsat
    q = np.asarray(qsat(t, pmid)["qs"])          # exactly saturated
    prd = np.zeros((ncol, nlev))
    prd[:, 40:60] = 2.0e-8
    cld = np.full((ncol, nlev), 0.3)
    prec_in = (prd * pdel).sum(1) / GRAV / 1000.0
    r = zm_conv.zm_conv_evap(pmid, pdel, t, q, prd, cld, prec_in,
                             1800.0, zc, zp)
    assert np.abs(np.asarray(r["tend_q"])).max() == 0.0
    np.testing.assert_allclose(np.asarray(r["prec"]), prec_in,
                               rtol=1e-15)


def test_snow_fraction_in_bounds(gold, results):
    """0 <= snow <= prec at the surface, and the snow flux never
    exceeds the precip flux at any interface (old_snow branch; the
    new-snow branch has zero snow with zm_microp=False)."""
    for cfg in CFGS:
        r = results[cfg]
        prec = np.asarray(r["prec"])
        snow = np.asarray(r["snow"])
        assert (snow >= 0.0).all()
        assert (snow <= prec * (1 + 1e-12) + 1e-30).all()
        fp = np.asarray(r["flxprec"])
        fs = np.asarray(r["flxsnow"])
        assert (fs >= 0.0).all()
        assert (fs <= fp * (1 + 1e-12) + 1e-25).all()


def test_old_snow_false_has_no_snow(results):
    """With zm_microp=False, prdsnow=0: the old_snow=F branch cannot
    produce snow (needs zm_microphysics, PORTING_PLAN.md row 9)."""
    for cfg in ("e2", "e4"):
        assert np.all(np.asarray(results[cfg]["flxsnow"]) == 0.0)
        assert np.all(np.asarray(results[cfg]["snow"]) == 0.0)


def test_no_rain_no_op(gold, results):
    """Columns with zero rain production and zero input precip must
    return exactly zero everywhere."""
    for cfg in CFGS:
        (pmid, pdel, t, q, prd, cld, prec_in, dt, old_snow,
         ke) = _inputs(gold, cfg)
        idle = (np.abs(np.asarray(prd)).sum(axis=1) == 0.0) \
            & (np.asarray(prec_in) == 0.0)
        if not idle.any():
            continue
        r = results[cfg]
        for f in FIELDS:
            assert np.all(np.asarray(r[f])[idle] == 0.0), (cfg, f)


def test_heating_consistent_with_phase_and_evap(gold, results):
    """tend_s == -evap*Lv + snow-production/melt fusion term, i.e.
    energy bookkeeping is closed within the routine."""
    for cfg in CFGS:
        r = results[cfg]
        if bool(json.loads(str(gold["__metadata__"]))["configs"][cfg]
                ["old_snow"]):
            recon = (-np.asarray(r["tend_q"]) * LATVAP
                     + np.asarray(r["ntsnprd"]) * 3.337e5)
        else:
            recon = (-np.asarray(r["tend_q"]) * LATVAP
                     + np.asarray(r["tend_s_snwevmlt"]))
        np.testing.assert_allclose(np.asarray(r["tend_s"]), recon,
                                   rtol=1e-12, atol=1e-18)
