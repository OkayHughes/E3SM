"""Tier-1 golden replay + Tier-0 properties for CLUBB slice D — the
Lscale/tau infrastructure (eam_jax.clubb.mixing_length, .stability,
.surface_varnce, .lscale_tau).

Golden: golden/clubb_lscale.npz — the REAL public routines
(compute_mixing_length, calc_brunt_vaisala_freq_sqd,
calc_stability_correction, term_wp2_splat/term_wp3_splat,
calc_surface_varnce) driven directly, plus drv_lscale_tau_segment, a
verbatim transcription of advance_clubb_core's inline
em->thvm->Lscale->tau->Kh->splat->surface->tau_N2 segment that the
generator validated BITWISE against the real advance_clubb_core
through its khzm/khzt diagnostics on all 25 cases.

MEASURED agreement (documented):
  BITWISE (replayed with eager op-by-op dispatch — under jax.jit,
  XLA's cross-op fusion introduces 1-ulp FMA-like contraction noise,
  measured <= 4.5e-16 rel on the same kernels): seg em,
  term_wp2/wp3_splat (incl. the exact (-0.0) EAMv3 C_wp2_splat=0
  tendencies), calc_surface_varnce direct (all six moments, all 44
  sweep cases), bv dry_t0 and moist variants.
  1-ulp (max rel <= 4e-16): thvm, sqrt_em_zt, bv dry_thvm, the
  segment surface moments (tau_zm(1) ulps through the splat guard).
  Lscale family: compute_mixing_length direct <= 1.0e-12 rel
  (83-96% of points bitwise); through the segment (thvm/sqrt_em 1-ulp
  inputs) lscale/tau/Kh/tau_N2 <= 2.8e-12 rel, stability_correction
  <= 4.7e-13.  The residue is the container-vs-XLA libm exp() in
  exp(-mu*dzm) feeding the sequential parcel recursion (the
  TKE-exhaustion stopping point is continuous across level
  boundaries, so branch flips do not amplify).  Asserted at
  rtol=2e-11 (7x measured), atol=0 — no field needed an atol floor.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import jax  # noqa: E402
import jax.numpy as jnp  # noqa: E402

from eam_jax.clubb import grid as cgrid  # noqa: E402
from eam_jax.clubb import lscale_tau as lt  # noqa: E402
from eam_jax.clubb import mixing_length as ml  # noqa: E402
from eam_jax.clubb import stability as stab  # noqa: E402
from eam_jax.clubb import surface_varnce as sv  # noqa: E402
from eam_jax.clubb.constants import (GRAV, RT_TOL, THL_TOL,  # noqa: E402
                                     W_TOL_SQD)

GOLDEN = ROOT / "golden"

SEG_IN = ["thlm", "rtm", "rcm", "wp2", "wp3", "up2", "vp2", "um",
          "vm", "p", "exner", "thv_ds_zt"]
SV_IN = ["upwp", "vpwp", "wpthlp", "wprtp", "um", "vm", "lup",
         "splat", "tau"]

RTOL_LSCALE = 2e-11     # measured <= 2.8e-12 (see module docstring)
RTOL_ULP = 1e-12        # measured <= 4e-16


@pytest.fixture(scope="module")
def g():
    path = GOLDEN / "clubb_lscale.npz"
    if not path.exists():
        pytest.skip("golden archive missing (run harness in container)")
    return np.load(path)


@pytest.fixture(scope="module")
def meta(g):
    return json.loads(str(g["meta"]))


@pytest.fixture(scope="module")
def tunables(g, meta):
    params = g["params"]
    lidx = meta["lscale_param_indices"]
    return dict(
        mu=meta["lscale_config"]["mu"],
        c_K=float(params[lidx["c_K"]]),
        taumax=float(params[lidx["taumax"]]),
        C_wp2_splat=float(params[lidx["C_wp2_splat"]]),
        lmin=float(g["lmin"]),
        lambda0_stability_coef=float(
            params[lidx["lambda0_stability_coef"]]),
        up2_vp2_factor=float(params[lidx["up2_vp2_factor"]]),
        T0=float(g["T0"]))


@pytest.fixture(scope="module")
def grid_f(g):
    return cgrid.setup_grid(g["zi"], g["zt"])


def test_config_is_eamv3(g, meta, tunables):
    """The golden was generated under the exact EAMv3 configuration
    this port encodes (flags read back from the Fortran modules)."""
    # drv_lscale_config slots (see clubb_driver.F90)
    assert list(g["lscale_flags"]) == [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
    cfg = meta["lscale_config"]
    assert cfg["l_stability_correct_tau_zm"] is True
    assert cfg["l_avg_Lscale"] is False          # perturbed avg DEAD
    assert cfg["l_use_C7_Richardson"] is False   # Cx_fnc dead
    assert cfg["l_use_C11_Richardson"] is False
    assert cfg["l_use_wp3_pr3"] is False
    assert tunables["mu"] == 0.0005 and tunables["c_K"] == 0.2
    assert tunables["taumax"] == 3600.0
    assert tunables["C_wp2_splat"] == 0.0
    assert tunables["lmin"] == 4.0 and tunables["T0"] == 300.0


# ------------------------------------------------- Tier-1: replays

@pytest.fixture(scope="module")
def ml_replay(g, grid_f):
    fn = jax.jit(ml.compute_mixing_length)
    outs = {k: [] for k in ["lscale", "lscale_up", "lscale_down"]}
    for c in range(g["ml_in_thvm"].shape[0]):
        ls, lu, ld = fn(grid_f, g["ml_in_thvm"][c], g["ml_in_thlm"][c],
                        g["ml_in_rtm"][c], g["ml_in_em"][c],
                        g["ml_lscale_max"][c], g["ml_in_p"][c],
                        g["ml_in_exner"][c], g["ml_in_thv_ds"][c],
                        g["ml_mu"][c], float(g["lmin"]))
        outs["lscale"].append(ls)
        outs["lscale_up"].append(lu)
        outs["lscale_down"].append(ld)
    return {k: np.stack(v) for k, v in outs.items()}


@pytest.mark.parametrize("field", ["lscale", "lscale_up",
                                   "lscale_down"])
def test_mixing_length_golden(g, ml_replay, field):
    np.testing.assert_allclose(ml_replay[field], g[f"ml_out_{field}"],
                               rtol=RTOL_LSCALE, atol=0.0)


def test_mixing_length_majority_bitwise(g, ml_replay):
    """Most points must stay bitwise (guards silent formula drift
    that would still pass rtol)."""
    for field in ["lscale", "lscale_up", "lscale_down"]:
        frac = (ml_replay[field] == g[f"ml_out_{field}"]).mean()
        assert frac > 0.7, (field, frac)


def _bv_fn(grid_f, t0, **kw):
    # eager (no jit): XLA jit fusion introduces 1-ulp (FMA-like)
    # contraction differences; op-by-op dispatch is BITWISE against
    # the -ffp-contract=off Fortran for these pure-arithmetic kernels
    return jax.vmap(
        lambda a, b, c, d, e, f: stab.calc_brunt_vaisala_freq_sqd(
            grid_f, a, b, c, d, e, f, t0, **kw))


@pytest.mark.parametrize("variant,kw", [
    ("dry_t0", {}),
    ("dry_thvm", {"l_use_thvm_in_bv_freq": True}),
    ("moist", {"l_brunt_vaisala_freq_moist": True}),
])
def test_brunt_vaisala_golden(g, grid_f, variant, kw):
    args = [jnp.asarray(g[f"bv_in_{k}"]) for k in
            ["thlm", "exner", "rtm", "rcm", "p", "thvm"]]
    got = np.asarray(_bv_fn(grid_f, float(g["T0"]), **kw)(*args))
    np.testing.assert_allclose(got, g[f"bv_out_{variant}"],
                               rtol=RTOL_ULP, atol=0.0)
    if variant in ("dry_t0", "moist"):    # measured bitwise (eager)
        assert (got == g[f"bv_out_{variant}"]).all()


def test_term_splat_golden_bitwise(g, grid_f):
    dt = float(g["dt"])
    for c in range(g["sp_in_wp2"].shape[0]):
        cs = float(g["sp_c_wp2_splat"][c])
        w2s = np.asarray(stab.term_wp2_splat(
            grid_f, cs, dt, g["sp_in_wp2"][c], g["sp_in_wp2_zt"][c],
            g["sp_in_tau_zm"][c]))
        w3s = np.asarray(stab.term_wp3_splat(
            grid_f, cs, dt, g["sp_in_wp2"][c], g["sp_in_wp3"][c],
            g["sp_in_tau_zt"][c]))
        assert np.array_equal(w2s, g["sp_out_wp2_splat"][c]), c
        assert np.array_equal(w3s, g["sp_out_wp3_splat"][c]), c
    # EAMv3 C_wp2_splat = 0: tendencies are exactly zero
    assert (g["sp_out_wp2_splat"][g["sp_c_wp2_splat"] == 0.0]
            == 0.0).all()
    # the nonzero sweep engaged the five/dt clip somewhere
    clip = 5.0 / dt
    nz_cases = g["sp_c_wp2_splat"] != 0.0
    assert (g["sp_out_wp2_splat"][nz_cases]
            == -g["sp_in_wp2"][nz_cases] * clip).any()


@pytest.fixture(scope="module")
def sv_replay(g, tunables):
    # eager vmap (no jit): bitwise; see _bv_fn note
    fn = jax.vmap(lambda *a: sv.calc_surface_varnce(
        *a, tunables["T0"], tunables["up2_vp2_factor"]))
    got = fn(*[jnp.asarray(g[f"sv_in_{k}"]) for k in SV_IN])
    return np.stack([np.asarray(x) for x in got], axis=1)


def test_surface_varnce_golden_bitwise(g, sv_replay):
    assert np.array_equal(sv_replay, g["sv_outs"])


@pytest.fixture(scope="module")
def seg_replay(g, grid_f, tunables):
    fn = jax.jit(lambda *a: lt.lscale_tau_segment(grid_f, *a,
                                                  tunables))
    dt = float(g["dt"])
    hd = float(g["host_dxy"])
    nseg, nz = g["seg_in_thlm"].shape
    seg = np.zeros((nseg, nz, 14))
    sfc = np.zeros((nseg, 6))
    extra = {"wp2_zt": np.zeros((nseg, nz)),
             "cx": np.zeros((nseg, nz))}
    for c in range(nseg):
        fl = g["seg_fluxes"][c]
        r = fn(dt, float(g["seg_sfc_elevation"][c]), hd, hd,
               fl[0], fl[1], fl[2], fl[3],
               *[g[f"seg_in_{k}"][c] for k in SEG_IN])
        for s in range(14):
            seg[c, :, s] = np.asarray(r[s])
        extra["wp2_zt"][c] = np.asarray(r.wp2_zt)
        extra["cx"][c] = np.asarray(r.cx_fnc_richardson)
        sfc[c] = [np.asarray(x) for x in
                  (r.wp2_sfc, r.up2_sfc, r.vp2_sfc, r.thlp2_sfc,
                   r.rtp2_sfc, r.rtpthlp_sfc)]
    return seg, sfc, extra


SEG_ULP_SLOTS = ["em", "thvm", "sqrt_em_zt"]
SEG_BITWISE_SLOTS = ["em", "wp2_splat", "wp3_splat"]


def test_segment_golden(g, meta, seg_replay):
    seg, _sfc, _extra = seg_replay
    for s, name in enumerate(meta["seg_slots"]):
        rtol = RTOL_ULP if name in SEG_ULP_SLOTS else RTOL_LSCALE
        np.testing.assert_allclose(seg[:, :, s],
                                   g["seg_outs"][:, :, s],
                                   rtol=rtol, atol=0.0, err_msg=name)
        if name in SEG_BITWISE_SLOTS:
            assert (seg[:, :, s] == g["seg_outs"][:, :, s]).all(), name


def test_segment_sfc_golden(g, meta, seg_replay):
    _seg, sfc, _extra = seg_replay
    for s, name in enumerate(meta["sfc_slots"]):
        np.testing.assert_allclose(sfc[:, s], g["seg_sfc_outs"][:, s],
                                   rtol=RTOL_ULP, atol=0.0,
                                   err_msg=name)
    # the elevated-surface case (last) takes the tolerance branch
    np.testing.assert_array_equal(
        sfc[-1], [W_TOL_SQD, W_TOL_SQD, W_TOL_SQD, THL_TOL ** 2,
                  RT_TOL ** 2, 0.0])


# ------------------------------------------------------ Tier-0

def test_tier0_lscale_bounds(g, ml_replay):
    """Lscale positive and capped by Lscale_max; up/down floored by
    the surface-layer lminh above the ghost level."""
    lminh = np.maximum(0.0, 500.0 - (np.asarray(g["zt"])
                                     - np.asarray(g["zi"])[0])) \
        * float(g["lmin"]) / 500.0
    for k in ["lscale", "lscale_up", "lscale_down"]:
        assert (ml_replay[k] >= 0.0).all()
    ls = ml_replay["lscale"]
    assert (ls[:, 1:] > 0.0).all()
    assert (ls <= g["ml_lscale_max"][:, None]).all()
    assert (ls == g["ml_lscale_max"][:, None]).any()   # cap exercised
    for k in ["lscale_up", "lscale_down"]:
        assert (ml_replay[k][:, 1:] >= lminh[None, 1:]).all()
        assert (ml_replay[k][:, 0] == 0.0).all()       # ghost level
    # geometric mean identity on interior levels
    up, dn = ml_replay["lscale_up"], ml_replay["lscale_down"]
    interior = np.minimum(np.sqrt(up * dn)[:, 1:-1],
                          g["ml_lscale_max"][:, None])
    np.testing.assert_allclose(ls[:, 1:-1], interior, rtol=1e-15)


def test_tier0_tau_positive_bounded(g, meta, seg_replay, tunables):
    seg, _sfc, extra = seg_replay
    slots = {n: i for i, n in enumerate(meta["seg_slots"])}
    tau_zt = seg[:, :, slots["tau_zt"]]
    tau_zm = seg[:, :, slots["tau_zm"]]
    tau_n2 = seg[:, :, slots["tau_n2_zm"]]
    sc = seg[:, :, slots["stability_correction"]]
    assert (tau_zt > 0.0).all() and (tau_zt <= tunables["taumax"]).all()
    assert (tau_zm > 0.0).all() and (tau_zm <= tunables["taumax"]).all()
    # stability correction: 1 <= sc <= 4, and tau_N2 = tau_zm/sc
    assert (sc >= 1.0).all() and (sc <= 4.0).all()
    assert (tau_n2 <= tau_zm * (1 + 1e-15)).all()
    np.testing.assert_array_equal(tau_n2, tau_zm / sc)
    # Kh composition identities (exact recomposition)
    np.testing.assert_array_equal(
        seg[:, :, slots["kh_zt"]],
        tunables["c_K"] * seg[:, :, slots["lscale"]]
        * seg[:, :, slots["sqrt_em_zt"]])
    # Cx_fnc_Richardson is dead in EAMv3
    assert (extra["cx"] == 0.0).all()
    assert (extra["wp2_zt"] >= W_TOL_SQD).all()


def test_tier0_brunt_vaisala_sign_matches_stability(grid_f):
    """Dry N^2 = (g/T0)*ddzt(thlm): stably stratified (thlm increasing
    with height) gives N^2 > 0 everywhere, unstable gives N^2 < 0, and
    uniform thlm gives exactly 0."""
    nz = 73
    z = np.linspace(0.0, 20000.0, nz)
    dummy = jnp.zeros(nz)
    for slope, sign in [(0.004, 1.0), (-0.004, -1.0)]:
        thlm = jnp.asarray(290.0 + slope * z)
        bv = np.asarray(stab.calc_brunt_vaisala_freq_sqd(
            grid_f, thlm, dummy, dummy, dummy, dummy, dummy, 300.0))
        assert (np.sign(bv) == sign).all()
    bv0 = np.asarray(stab.calc_brunt_vaisala_freq_sqd(
        grid_f, jnp.full(nz, 300.0), dummy, dummy, dummy, dummy,
        dummy, 300.0))
    assert (bv0 == 0.0).all()


def test_tier0_surface_variance_thresholds(g, sv_replay):
    """Floors/realizability of the surface moments across the flux
    sweep."""
    wp2, up2, vp2 = sv_replay[:, 0], sv_replay[:, 1], sv_replay[:, 2]
    thlp2, rtp2, rtpthlp = (sv_replay[:, 3], sv_replay[:, 4],
                            sv_replay[:, 5])
    assert (wp2 >= W_TOL_SQD).all()
    assert (thlp2 >= THL_TOL ** 2).all()
    assert (rtp2 >= RT_TOL ** 2).all()
    assert (up2 > 0.0).all() and (vp2 > 0.0).all()
    # rt/thl covariance realizability: |rtpthlp| = 0.2a|x||y| while
    # rtp2*thlp2 >= (0.4a x^2)(0.4a y^2) => |corr| <= 0.5
    assert (np.abs(rtpthlp)
            <= np.sqrt(rtp2 * thlp2) * (1 + 1e-12)).all()
    # zero-flux cases give exactly the floors
    zero = (g["sv_in_wpthlp"] == 0.0) & (g["sv_in_wprtp"] == 0.0)
    assert (thlp2[zero] == THL_TOL ** 2).all()
    assert (rtp2[zero] == RT_TOL ** 2).all()
    assert (rtpthlp[zero] == 0.0).all()


def test_tier0_lscale_monotone_in_stability(grid_f, g, tunables):
    """The parcel argument: with identical TKE, a more stably
    stratified column exhausts ASCENDING parcels sooner, so Lscale_up
    can only shrink as d(thv)/dz increases (level-by-level).  The
    combined Lscale = sqrt(up*down) is NOT asserted: the descending
    parcels' nonlocal min-altitude smoothing makes Lscale_down (and
    hence the geometric mean) non-monotone at isolated levels."""
    zt = np.asarray(g["zt"])
    zi = np.asarray(g["zi"])
    nz = zt.size
    z = np.maximum(zt, 0.0)
    p = 1.0e5 * np.exp(-z / 8000.0)
    p[0] = p[1]
    exner = (p / 1.0e5) ** (287.042 / 1004.64)
    em = jnp.asarray(0.3 * np.exp(-np.maximum(zi, 0.0) / 1500.0)
                     + 1e-3)
    rtm = jnp.full(nz, 2e-3)
    prev_up = None
    for slope in [0.001, 0.002, 0.004, 0.008]:
        thlm = jnp.asarray(290.0 + slope * z)
        thv_ds = thlm
        thvm = thlm + 0.61 * thv_ds * rtm
        _ls, lu, _ld = ml.compute_mixing_length(
            grid_f, thvm, thlm, rtm, em, 25000.0, jnp.asarray(p),
            jnp.asarray(exner), thv_ds, tunables["mu"],
            tunables["lmin"])
        lu = np.asarray(lu)
        if prev_up is not None:
            assert (lu <= prev_up * (1 + 1e-12)).all()
        prev_up = lu


def test_tier0_lscale_monotone_in_tke(grid_f, g, tunables):
    """More initial parcel TKE cannot shorten the ascent: Lscale_up
    and Lscale_down are nondecreasing in em."""
    zt = np.asarray(g["zt"])
    nz = zt.size
    z = np.maximum(zt, 0.0)
    p = 1.0e5 * np.exp(-z / 8000.0)
    p[0] = p[1]
    exner = (p / 1.0e5) ** (287.042 / 1004.64)
    thlm = jnp.asarray(288.0 + 0.003 * z)
    rtm = jnp.full(nz, 3e-3)
    thvm = thlm + 0.61 * thlm * rtm
    prev = None
    for em0 in [0.01, 0.05, 0.2, 1.0]:
        ls, lu, ld = ml.compute_mixing_length(
            grid_f, thvm, thlm, rtm, jnp.full(nz, em0), 25000.0,
            jnp.asarray(p), jnp.asarray(exner), thlm, tunables["mu"],
            tunables["lmin"])
        cur = (np.asarray(lu), np.asarray(ld))
        if prev is not None:
            assert (cur[0] >= prev[0] * (1 - 1e-12)).all()
            assert (cur[1] >= prev[1] * (1 - 1e-12)).all()
        prev = cur


def test_tier0_splat_properties(grid_f, g):
    """wp2_splat <= 0 (a sink), exactly zero for C=0, and wp3_splat
    opposes wp3."""
    dt = float(g["dt"])
    nz = g["zi"].size
    rng = np.random.default_rng(5)
    wp2 = jnp.asarray(0.5 * rng.random(nz) + W_TOL_SQD)
    wp2_zt = jnp.maximum(cgrid.zm2zt(grid_f, wp2), W_TOL_SQD)
    wp3 = jnp.asarray(rng.standard_normal(nz) * 0.1)
    tau = jnp.full(nz, 900.0)
    w2s = np.asarray(stab.term_wp2_splat(grid_f, 2.0, dt, wp2,
                                         wp2_zt, tau))
    w3s = np.asarray(stab.term_wp3_splat(grid_f, 2.0, dt, wp2, wp3,
                                         tau))
    assert (w2s <= 0.0).all()
    assert (w3s * np.asarray(wp3) <= 0.0).all()
    assert (np.asarray(stab.term_wp2_splat(grid_f, 0.0, dt, wp2,
                                           wp2_zt, tau)) == 0.0).all()
    # clip: |wp2_splat| <= (5/dt)*wp2
    assert (np.abs(w2s) <= 5.0 / dt * np.asarray(wp2)
            * (1 + 1e-15)).all()


def test_tier0_stability_correction_neutral_is_one(grid_f):
    """Uniform thlm (N^2 = 0 exactly) gives stability_correction == 1
    everywhere, so tau_N2_zm == tau_zm."""
    nz = 73
    thlm = jnp.full(nz, 300.0)
    ones = jnp.ones(nz)
    sc = np.asarray(stab.calc_stability_correction(
        grid_f, thlm, 100.0 * ones, 0.3 * ones, ones, 1e-3 * ones,
        jnp.zeros(nz), 9e4 * ones, thlm, 0.03, 300.0))
    assert (sc == 1.0).all()
