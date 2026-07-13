"""Tier-1 golden replay + Tier-0 properties for
eam_jax.clubb.advance_xp2_xpyp (CLUBB slice C: the prognostic advance
of rt'2, thl'2, rt'thl', u'2, v'2 plus clip_covars_denom).

Golden: golden/clubb_xp2.npz — the REAL public advance_xp2_xpyp
(eam_clubb_f, the full unmodified CLUBB stack) driven directly over 40
cases: 16 whose prognostic moments/Kh are one-step
advance_clubb_core-advanced states of the slice-B regime columns, and
24 synthetic stress columns (degenerate variances, Cauchy-Schwarz
violations, sign-alternating wp3_on_wp2 for both upwind arms, the
0.5*rtm^2 rtp2 cap, the 1000 m2/s2 up2/vp2 cap, hole-filling baited by
strongly negative forcings, advection/dissipation-heavy, randomized);
plus clip_covars_denom on the post-xp2 states.

MEASURED agreement (documented; asserted at rtol=1e-12, atol=0):
  clip_covars_denom: BITWISE on all four covariances, all 40 cases.
  advance_xp2_xpyp: 59-73% of points bitwise; the rest a few ulp of
  their own value — max rel: rtp2 3.5e-15, thlp2 2.4e-15, up2 1.2e-15,
  vp2 1.1e-15, rtpthlp 8.5e-14 (at ~4e-8 leading-edge magnitudes where
  cancelling turbulent-production terms lose common bits; abs
  6.5e-19).  The residue is XLA FMA contraction inside the LHS/RHS
  term assemblies (the Fortran objects are -ffp-contract=off)
  propagated through the dgtsv solve; no field needed a loosened
  (atol) floor.
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

from eam_jax.clubb import advance_xp2_xpyp as axp  # noqa: E402
from eam_jax.clubb import grid as cgrid  # noqa: E402
from eam_jax.clubb.constants import (MAX_MAG_CORRELATION,  # noqa: E402
                                     RT_TOL, THL_TOL, W_TOL_SQD)

GOLDEN = ROOT / "golden"

# advance_xp2_xpyp port inputs, in signature order (golden keys carry
# the in_ prefix; the golden additionally records the DEAD Fortran
# arguments skw_zm/cloud_frac/wprtp2/wpthlp2/wprtpthlp/wp3/lscale/
# wp2_zt/wp3_on_wp2_zt/rho_ds_zt, which the port does not take)
XP2_IN = ["tau_zm", "wm_zm", "rtm", "wprtp", "thlm", "wpthlp",
          "wpthvp", "um", "vm", "wp2", "upwp", "vpwp", "sigma_sqd_w",
          "kh_zt", "rtp2_forcing", "thlp2_forcing", "rtpthlp_forcing",
          "rho_ds_zm", "invrs_rho_ds_zm", "thv_ds_zm", "wp3_on_wp2",
          "wp2_splat", "rtp2", "thlp2", "rtpthlp", "up2", "vp2"]

XP2_OUT = ["rtp2", "thlp2", "rtpthlp", "up2", "vp2"]

CC_IN = ["wp2", "rtp2", "thlp2", "up2", "vp2", "wprtp", "wpthlp",
         "upwp", "vpwp"]
CC_OUT = ["wprtp", "wpthlp", "upwp", "vpwp"]


@pytest.fixture(scope="module")
def gxp2():
    path = GOLDEN / "clubb_xp2.npz"
    if not path.exists():
        pytest.skip("golden archive missing (run harness in container)")
    return np.load(path)


@pytest.fixture(scope="module")
def tunables(gxp2):
    meta = json.loads(str(gxp2["meta"]))
    params = gxp2["params"]
    idx = meta["param_indices"]
    xidx = meta["xp2_param_indices"]
    return dict(
        C2rt=float(params[idx["C2rt"]]),
        C2thl=float(params[idx["C2thl"]]),
        C2rtthl=float(params[idx["C2rtthl"]]),
        beta=float(params[idx["beta"]]),
        C4=float(params[xidx["C4"]]),
        C5=float(params[xidx["C5"]]),
        C14=float(params[xidx["C14"]]),
        c_K2=float(params[xidx["c_K2"]]),
        c_K9=float(params[xidx["c_K9"]]),
    )


@pytest.fixture(scope="module")
def grid_f(gxp2):
    return cgrid.setup_grid(gxp2["zi"], gxp2["zt"])


def _advance_fn(gxp2, grid, tun, **kw):
    dt = float(gxp2["dt"])
    nu2 = jnp.asarray(gxp2["nu2_vert_res_dep"])
    nu9 = jnp.asarray(gxp2["nu9_vert_res_dep"])

    def one(*args):
        return axp.advance_xp2_xpyp(grid, dt, *args, nu2, nu9, tun,
                                    **kw)

    return jax.jit(jax.vmap(one))


@pytest.fixture(scope="module")
def replay(gxp2, grid_f, tunables):
    fn = _advance_fn(gxp2, grid_f, tunables)
    ins = [jnp.asarray(gxp2[f"in_{k}"]) for k in XP2_IN]
    outs = fn(*ins)
    return {k: np.asarray(v) for k, v in zip(XP2_OUT, outs)}


@pytest.fixture(scope="module")
def replay_nofill(gxp2, grid_f, tunables):
    fn = _advance_fn(gxp2, grid_f, tunables, l_hole_fill=False)
    ins = [jnp.asarray(gxp2[f"in_{k}"]) for k in XP2_IN]
    outs = fn(*ins)
    return {k: np.asarray(v) for k, v in zip(XP2_OUT, outs)}


# ------------------------------------------------------ Tier-1: golden

@pytest.mark.parametrize("field", XP2_OUT)
def test_xp2_golden(gxp2, replay, field):
    ref = gxp2[f"out_{field}"]
    got = replay[field]
    np.testing.assert_allclose(got, ref, rtol=1e-12, atol=0.0)


def test_xp2_agreement_documented(gxp2, replay):
    """Guard the measured agreement quality: a majority of points must
    stay bitwise (catches silent wholesale formula drift that still
    passes rtol)."""
    for field in XP2_OUT:
        frac = (replay[field] == gxp2[f"out_{field}"]).mean()
        assert frac > 0.5, (field, frac)


def test_clip_covars_denom_golden_bitwise(gxp2):
    fn = jax.jit(jax.vmap(axp.clip_covars_denom))
    outs = fn(*[jnp.asarray(gxp2[f"cc_in_{k}"]) for k in CC_IN])
    for name, got in zip(CC_OUT, outs):
        assert np.array_equal(np.asarray(got), gxp2[f"cc_out_{name}"]), \
            name


def test_golden_branch_coverage(gxp2):
    """The golden set must actually exercise every clip (mirrors the
    generator's coverage asserts, re-checked on the archive)."""
    rtp2 = gxp2["out_rtp2"]
    assert (rtp2 == RT_TOL ** 2).any()
    assert (gxp2["out_thlp2"] == THL_TOL ** 2).any()
    assert (gxp2["out_up2"] == 1000.0).any()
    assert (rtp2 == 0.5 * gxp2["in_rtm"] ** 2).any()
    bound = MAX_MAG_CORRELATION * np.sqrt(rtp2 * gxp2["out_thlp2"])
    assert (np.isclose(np.abs(gxp2["out_rtpthlp"]), bound, rtol=1e-14)
            & (bound > 0)).any()
    # clip_covars_denom clipped something
    assert any((gxp2[f"cc_out_{k}"] != gxp2[f"cc_in_{k}"]).any()
               for k in CC_OUT)


def test_hole_filling_fires_on_goldens(gxp2, replay, replay_nofill):
    """Branch coverage: pos_definite_variances must matter for at
    least one golden case (the hole-bait family), i.e. disabling
    l_hole_fill changes the replay."""
    changed = any((replay[k] != replay_nofill[k]).any() for k in XP2_OUT)
    assert changed


# ------------------------------------------------------ Tier-0

def test_tier0_realizability_post_clip(gxp2, replay):
    """Positivity/bounds after the clipping sequence, on the replayed
    outputs of every golden case."""
    nz = replay["rtp2"].shape[1]
    rtm = gxp2["in_rtm"]
    for c in range(replay["rtp2"].shape[0]):
        rtp2, thlp2 = replay["rtp2"][c], replay["thlp2"][c]
        up2, vp2 = replay["up2"][c], replay["vp2"][c]
        # clip_variance floors every level but the top; the large-rtp2
        # cap (all levels) may pull rtp2 below rt_tol^2 where rtm is
        # tiny
        floor_rt = np.minimum(RT_TOL ** 2, 0.5 * rtm[c] ** 2)
        assert (rtp2[:-1] >= floor_rt[:-1]).all()
        assert (rtp2 >= 0.0).all()
        assert (rtp2 <= 0.5 * rtm[c] ** 2 + 1e-30).all()
        assert (thlp2[:-1] >= THL_TOL ** 2).all()
        assert (up2[:-1] >= W_TOL_SQD).all() and (up2 <= 1000.0).all()
        assert (vp2[:-1] >= W_TOL_SQD).all() and (vp2 <= 1000.0).all()
        # covariance realizability on the clipped interior
        bound = MAX_MAG_CORRELATION * np.sqrt(rtp2 * thlp2)
        r = replay["rtpthlp"][c]
        assert (np.abs(r[1:nz - 1]) <= bound[1:nz - 1]
                * (1.0 + 1e-15)).all()


def _quiescent_inputs(nz, x_rt, x_thl, x_uv, wp2_val, tau=600.0,
                      thlm_grad=None, zt=None, wpthlp_val=0.0):
    """Inputs with every transport/production term off (wm=0, Kh=0,
    fluxes 0, wp3_on_wp2=0 -> zero ta) so only the dp1/pr1 dissipation
    balance remains; optionally a linear thlm(z) + constant wpthlp for
    a nonzero turbulent-production case."""
    z = np.zeros(nz)
    thlm = np.full(nz, 290.0) if thlm_grad is None \
        else 290.0 + thlm_grad * zt
    ins = dict(
        tau_zm=np.full(nz, tau), wm_zm=z, rtm=np.full(nz, 8e-3),
        wprtp=z, thlm=thlm, wpthlp=np.full(nz, wpthlp_val),
        wpthvp=z, um=np.full(nz, 5.0), vm=np.full(nz, -3.0),
        wp2=np.full(nz, wp2_val), upwp=z, vpwp=z,
        sigma_sqd_w=np.full(nz, 0.4), kh_zt=z, rtp2_forcing=z,
        thlp2_forcing=z, rtpthlp_forcing=z,
        rho_ds_zm=np.full(nz, 1.0), invrs_rho_ds_zm=np.full(nz, 1.0),
        thv_ds_zm=np.full(nz, 290.0), wp3_on_wp2=z, wp2_splat=z,
        rtp2=np.full(nz, x_rt), thlp2=np.full(nz, x_thl),
        rtpthlp=z, up2=np.full(nz, x_uv), vp2=np.full(nz, x_uv))
    return ins


def test_tier0_dissipation_threshold_fixed_point(gxp2, grid_f,
                                                 tunables):
    """With all transport and production off, the dissipation terms
    damp each variance toward its threshold — so a state sitting AT
    the thresholds (rtp2=rt_tol^2, thlp2=thl_tol^2,
    up2=vp2=wp2=w_tol^2) is a fixed point of the advance (the up2/vp2
    dp1+pr1 combination has the same fixed point:
    ((2C4+C14)/3) x* = (2/3)(C4-C14) w_tol^2 + C14 w_tol^2  =>
    x* = w_tol^2)."""
    nz = gxp2["zi"].size
    dt = float(gxp2["dt"])
    ins = _quiescent_inputs(nz, RT_TOL ** 2, THL_TOL ** 2, W_TOL_SQD,
                            W_TOL_SQD)
    nu0 = jnp.zeros(nz)
    outs = axp.advance_xp2_xpyp(
        grid_f, dt, *[jnp.asarray(ins[k]) for k in XP2_IN], nu0, nu0,
        tunables)
    for name, ref in zip(XP2_OUT,
                         [RT_TOL ** 2, THL_TOL ** 2, 0.0, W_TOL_SQD,
                          W_TOL_SQD]):
        np.testing.assert_allclose(np.asarray(outs[XP2_OUT.index(name)]),
                                   ref, rtol=1e-12, atol=1e-25,
                                   err_msg=name)


def test_tier0_production_dissipation_balance(gxp2, grid_f, tunables):
    """Steady state with production == dissipation: a linear thlm(z)
    and constant wpthlp give a uniform turbulent production
    P = -2*wpthlp*dthlm/dz; the steady thlp2 is
    thl_tol^2 + P*tau/C2thl, and feeding that state back through the
    advance must return it unchanged (interior levels; the top level
    is pinned to thl_tol^2 by the boundary condition)."""
    nz = gxp2["zi"].size
    dt = float(gxp2["dt"])
    tau = 600.0
    s = -0.005          # dthlm/dz [K/m]
    w0 = 0.02           # wpthlp [K m/s]  -> P = -2*w0*s > 0
    P = -2.0 * w0 * s
    x_star = THL_TOL ** 2 + P * tau / tunables["C2thl"]
    ins = _quiescent_inputs(nz, RT_TOL ** 2, x_star, W_TOL_SQD,
                            W_TOL_SQD, tau=tau, thlm_grad=s,
                            zt=np.asarray(gxp2["zt"]), wpthlp_val=w0)
    nu0 = jnp.zeros(nz)
    outs = axp.advance_xp2_xpyp(
        grid_f, dt, *[jnp.asarray(ins[k]) for k in XP2_IN], nu0, nu0,
        tunables)
    thlp2 = np.asarray(outs[1])
    np.testing.assert_allclose(thlp2[1:-1], x_star, rtol=1e-12)
    assert thlp2[-1] == THL_TOL ** 2   # top BC pins to threshold


def test_tier0_dead_path_exact_noop(gxp2, grid_f, tunables):
    """A comfortable state (variances well above every threshold, no
    correlation violations, no holes) must be BITWISE identical with
    and without hole filling, and clip_covar/clip_variance must be
    exact no-ops on in-bound data."""
    nz = gxp2["zi"].size
    dt = float(gxp2["dt"])
    rng = np.random.default_rng(7)
    ins = _quiescent_inputs(nz, 1e-7, 0.5, 0.3, 0.5)
    ins["wm_zm"] = 0.01 * rng.standard_normal(nz)
    ins["kh_zt"] = np.full(nz, 10.0)
    nu2 = jnp.asarray(gxp2["nu2_vert_res_dep"])
    nu9 = jnp.asarray(gxp2["nu9_vert_res_dep"])
    args = [jnp.asarray(ins[k]) for k in XP2_IN]
    with_fill = axp.advance_xp2_xpyp(grid_f, dt, *args, nu2, nu9,
                                     tunables)
    without = axp.advance_xp2_xpyp(grid_f, dt, *args, nu2, nu9,
                                   tunables, l_hole_fill=False)
    for name, a, b in zip(XP2_OUT, with_fill, without):
        assert np.array_equal(np.asarray(a), np.asarray(b)), name

    # clip helpers: exact no-op on in-bound data
    xp2 = jnp.asarray(1e-3 + 1e-4 * rng.random(nz))
    assert np.array_equal(np.asarray(axp.clip_variance(1e-6, xp2)),
                          np.asarray(xp2))
    xpyp = 0.5 * jnp.sqrt(xp2 * xp2)
    assert np.array_equal(
        np.asarray(axp.clip_covar(MAX_MAG_CORRELATION, xp2, xp2, xpyp)),
        np.asarray(xpyp))
    cc = axp.clip_covars_denom(xp2, xp2, xp2, xp2, xp2, xpyp, xpyp,
                               xpyp, xpyp)
    for got in cc:
        assert np.array_equal(np.asarray(got), np.asarray(xpyp))


def test_tier0_fill_holes_conserves_mass(gxp2, grid_f):
    """fill_holes_vertical ("zm") conserves the density-weighted
    integral over the hole-filling range (Fortran levels 2..nz-1) and
    never touches the surface or top levels."""
    nz = gxp2["zi"].size
    rng = np.random.default_rng(3)
    rho = np.asarray(gxp2["in_rho_ds_zm"][0])
    field = 1e-3 * (1.0 + rng.random(nz))
    field[10:14] = -2e-4        # a hole
    field[40] = -1e-5           # another
    tol = 1e-8
    filled = np.asarray(axp.fill_holes_vertical_zm(
        grid_f, 2, tol, jnp.asarray(rho), jnp.asarray(field)))
    assert filled[0] == field[0] and filled[-1] == field[-1]
    dzm = np.asarray(grid_f.dzm)
    m0 = np.sum(rho[1:nz - 1] * dzm[1:nz - 1] * field[1:nz - 1])
    m1 = np.sum(rho[1:nz - 1] * dzm[1:nz - 1] * filled[1:nz - 1])
    np.testing.assert_allclose(m1, m0, rtol=1e-12)
    # holes were actually reduced
    assert filled[10:14].min() > field[10:14].min()


def test_tier0_clip_covar_lands_on_bound():
    """Violating covariances are clipped exactly onto
    +/- max_mag_corr*sqrt(xp2*yp2) on interior levels."""
    nz = 12
    xp2 = jnp.full(nz, 4.0)
    yp2 = jnp.full(nz, 9.0)
    xpyp = jnp.asarray([(-1.0) ** k * 50.0 for k in range(nz)])
    got = np.asarray(axp.clip_covar(MAX_MAG_CORRELATION, xp2, yp2,
                                    xpyp))
    bound = MAX_MAG_CORRELATION * np.sqrt(4.0 * 9.0)
    assert got[0] == 50.0 and got[-1] == -50.0  # boundaries untouched
    np.testing.assert_array_equal(np.abs(got[1:-1]), bound)
    np.testing.assert_array_equal(np.sign(got[1:-1]),
                                  np.sign(np.asarray(xpyp)[1:-1]))
