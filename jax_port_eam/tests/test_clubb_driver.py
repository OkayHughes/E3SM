"""Tier-1 golden replay + Tier-0 properties for
eam_jax.clubb.pdf_closure_driver (CLUBB slice B: the zt+zm double
pdf_closure call, trapezoidal-rule vertical averaging, clip_rcm,
compute_cloud_cover, l_use_cloud_cover substitution).

Golden: golden/clubb_pdf_driver.npz — the verbatim-extracted Fortran
pdf_closure_driver over 40 regime-sweeping columns (inputs on native
grids: moments on zm, means on zt; incl. a one-level-dry-notch family
that triggers clip_rcm and the compute_cloud_cover partial-fill
branches) plus 8 "adv" cases whose inputs are real
advance_clubb_core-advanced states.  The extraction itself was
validated BITWISE against the public advance_clubb_core during
generation (see harness/gen_clubb_golden.py).

MEASURED agreement (both input families): everything is at the
1-2-ulp level of its own scale.  The only abs errors above ~1e-15
are (a) 1-ulp FMA-contraction differences in the linear grid
interpolations at thl's ~800 K scale (thlm_zm/thl_1/thl_2 ~1e-13
abs = 1.4e-16 rel; the Fortran objects are -ffp-contract=off, XLA
contracts), and (b) the slice-A erf/exp libm tails amplified by
knife-edge zeta sensitivity in the second (zm) pdf_closure call —
1-ulp input differences in chi/stdev_chi divided by tiny stdev move
zeta by ~1e-11, so cloud_frac/ice_supersat_frac leading edges differ
by up to ~1e-11..1e-12 abs (their values there are <=1e-6), and the
cloud-cover division cf/vert_cloud_frac carries that through.  All
atol floors below are ~10x the measured maxima; rtol stays 1e-12.
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
from eam_jax.clubb import pdf_closure as pc  # noqa: E402
from eam_jax.clubb import pdf_closure_driver as pcd  # noqa: E402
from eam_jax.clubb.constants import EPSILON_R8, RC_TOL  # noqa: E402
from tests.test_clubb import PDFP_SLOTS  # noqa: E402

GOLDEN = ROOT / "golden"

DRIVER_INPUTS = ["wprtp", "thlm", "wpthlp", "rtp2", "rtp3", "thlp2",
                 "thlp3", "rtpthlp", "wp2", "wp3", "wm_zm", "wm_zt",
                 "um", "up2", "upwp", "vm", "vp2", "vpwp", "p",
                 "exner", "thv_ds_zm", "thv_ds_zt", "rtm"]

OUTS_ATOL = {
    "rcm": 5e-17, "cloud_frac": 2e-12, "ice_supersat_frac": 3e-11,
    "wprcp": 1e-16, "sigma_sqd_w": 5e-16, "wpthvp": 1e-13,
    "wp2thvp": 1e-13, "rtpthvp": 2e-15, "thlpthvp": 1e-13,
    "rc_coef": 0.0, "rcm_in_layer": 5e-17, "cloud_cover": 2e-12,
    "rcp2_zt": 5e-21, "thlprcp": 1e-16, "rc_coef_zm": 0.0,
    "wp2rtp": 1e-17, "wp2thlp": 1e-15, "wp2rcp": 5e-18,
    "rtprcp": 1e-18, "rcp2": 5e-19, "uprcp": 1e-16, "vprcp": 1e-16,
    "cloud_frac_zm": 1e-11, "ice_supersat_frac_zm": 1e-10,
    "rtm_zm": 1e-18, "thlm_zm": 2e-12, "rcm_zm": 5e-16,
    "rcm_supersat_adj": 0.0, "sigma_sqd_w_zt": 5e-16,
}

PDFP_ATOL = {
    "w_1": 3e-15, "w_2": 2e-15, "varnce_w_1": 3e-16,
    "varnce_w_2": 3e-16, "rt_1": 2e-17, "rt_2": 2e-17,
    "varnce_rt_1": 6e-20, "varnce_rt_2": 6e-20, "thl_1": 2e-12,
    "thl_2": 2e-12, "varnce_thl_1": 2e-14, "varnce_thl_2": 2e-14,
    "corr_rt_thl_1": 2e-14, "corr_rt_thl_2": 2e-14,
    "alpha_thl": 1e-15, "alpha_rt": 1e-15, "crt_1": 2e-14,
    "crt_2": 2e-14, "cthl_1": 1e-17, "cthl_2": 1e-17, "chi_1": 3e-16,
    "chi_2": 3e-16, "stdev_chi_1": 1e-16, "stdev_chi_2": 1e-16,
    "stdev_eta_1": 1e-16, "stdev_eta_2": 1e-16,
    "covar_chi_eta_1": 1e-18, "covar_chi_eta_2": 1e-18,
    "corr_chi_eta_1": 2e-13, "corr_chi_eta_2": 2e-13,
    "rsatl_1": 1e-15, "rsatl_2": 1e-15, "rc_1": 3e-16, "rc_2": 3e-16,
    "cloud_frac_1": 1e-11, "cloud_frac_2": 1e-11, "mixt_frac": 2e-15,
    "ice_supersat_frac_1": 2e-10, "ice_supersat_frac_2": 2e-10,
}


@pytest.fixture(scope="module")
def gdrv():
    path = GOLDEN / "clubb_pdf_driver.npz"
    if not path.exists():
        pytest.skip("golden archive missing (run harness in container)")
    return np.load(path)


@pytest.fixture(scope="module")
def tunables(gdrv):
    meta = json.loads(str(gdrv["meta"]))
    params = gdrv["params"]
    idx = meta["param_indices"]
    return dict(
        beta=float(params[idx["beta"]]),
        mf_max_mag=pc.mixt_frac_max_mag(
            float(params[idx["Skw_max_mag"]])),
        gamma_coef=float(params[idx["gamma_coef"]]),
        gamma_coefb=float(params[idx["gamma_coefb"]]),
        gamma_coefc=float(params[idx["gamma_coefc"]]),
    ), meta["outs_slots"]


@pytest.fixture(scope="module")
def grid_f(gdrv):
    return cgrid.setup_grid(gdrv["zi"], gdrv["zt"])


def _replay(gdrv, grid, tun, prefix):
    ins = [jnp.asarray(gdrv[f"{prefix}{k}"]) for k in DRIVER_INPUTS]
    fn = jax.vmap(lambda *a: pcd.pdf_closure_driver(grid, *a, **tun),
                  in_axes=0)
    outs, pz, pm = fn(*ins)
    return ({k: np.asarray(v) for k, v in outs.items()},
            {k: np.asarray(v) for k, v in pz.items()},
            {k: np.asarray(v) for k, v in pm.items()})


@pytest.fixture(scope="module")
def replay_syn(gdrv, grid_f, tunables):
    return _replay(gdrv, grid_f, tunables[0], "in_")


@pytest.fixture(scope="module")
def replay_adv(gdrv, grid_f, tunables):
    return _replay(gdrv, grid_f, tunables[0], "adv_in_")


# ------------------------------------------------------ Tier-1: golden
@pytest.mark.parametrize("family", ["syn", "adv"])
def test_driver_outs_golden(gdrv, tunables, replay_syn, replay_adv,
                            family):
    slots = tunables[1]
    outs = (replay_syn if family == "syn" else replay_adv)[0]
    ref = gdrv["outs"] if family == "syn" else gdrv["adv_outs"]
    for j, name in enumerate(slots):
        np.testing.assert_allclose(
            outs[name], ref[:, :, j], rtol=1e-12,
            atol=OUTS_ATOL[name], err_msg=f"{family}:{name}")


@pytest.mark.parametrize("family", ["syn", "adv"])
@pytest.mark.parametrize("which", ["zt", "zm"])
def test_driver_pdf_params_golden(gdrv, replay_syn, replay_adv,
                                  family, which):
    rep = replay_syn if family == "syn" else replay_adv
    pdfp = rep[1] if which == "zt" else rep[2]
    key = ("" if family == "syn" else "adv_") + f"pdfp_{which}"
    ref = gdrv[key]
    for j, name in enumerate(PDFP_SLOTS):
        np.testing.assert_allclose(
            pdfp[name], ref[:, :, j], rtol=1e-12,
            atol=PDFP_ATOL.get(name, 0.0),
            err_msg=f"{family}:{which}:{name}")


def test_rtm_passes_through(gdrv):
    # l_rtm_nudge = .false.: the Fortran driver's intent(inout) rtm is
    # bitwise unchanged (which is why the port does not return it).
    np.testing.assert_array_equal(gdrv["rtm_out"], gdrv["in_rtm"])
    np.testing.assert_array_equal(gdrv["adv_rtm_out"],
                                  gdrv["adv_in_rtm"])


def test_clip_rcm_coverage(gdrv, grid_f, tunables):
    """The dry-notch family (fam 4) must exercise the clip_rcm branch
    (rtm < rcm after the trapezoidal rule) at least once, so the
    golden actually validates that path."""
    tun = tunables[0]
    hits = []
    orig = pcd.clip_rcm

    def counting(rtm, rcm):
        hits.append(np.asarray(rtm < rcm))
        return orig(rtm, rcm)

    pcd.clip_rcm = counting
    try:
        for c in range(gdrv["in_rtm"].shape[0]):
            args = [jnp.asarray(gdrv[f"in_{k}"][c])
                    for k in DRIVER_INPUTS]
            pcd.pdf_closure_driver(grid_f, *args, **tun)
    finally:
        pcd.clip_rcm = orig
    assert np.array(hits).sum() >= 1


# ----------------------------------------------- Tier-0: trapezoid rule
def test_tier0_trapezoid_properties(grid_f):
    gr = grid_f
    nz = gr.zt.shape[0]
    rng = np.random.default_rng(11)
    f_zt = jnp.asarray(rng.normal(size=nz))
    f_zm = jnp.asarray(rng.normal(size=nz))

    # constant in == constant out (the weights sum to exactly 1)
    const = jnp.full(nz, 7.25)
    np.testing.assert_allclose(
        np.asarray(pcd.trapezoid_zt(gr, const, const)), 7.25,
        rtol=1e-15)
    np.testing.assert_allclose(
        np.asarray(pcd.trapezoid_zm(gr, const, const)), 7.25,
        rtol=1e-15)

    # each interior level is a convex combination of the three
    # surrounding values (positive weights)
    tz = np.asarray(pcd.trapezoid_zt(gr, f_zt, f_zm))
    lo = np.minimum(np.asarray(f_zt)[1:],
                    np.minimum(np.asarray(f_zm)[1:],
                               np.asarray(f_zm)[:-1]))
    hi = np.maximum(np.asarray(f_zt)[1:],
                    np.maximum(np.asarray(f_zm)[1:],
                               np.asarray(f_zm)[:-1]))
    assert (tz[1:] >= lo - 1e-12).all() and (tz[1:] <= hi + 1e-12).all()
    assert tz[0] == np.asarray(f_zt)[0]  # boundary copies zt value

    tm = np.asarray(pcd.trapezoid_zm(gr, f_zm, f_zt))
    lo = np.minimum(np.asarray(f_zm)[1:-1],
                    np.minimum(np.asarray(f_zt)[1:-1],
                               np.asarray(f_zt)[2:]))
    hi = np.maximum(np.asarray(f_zm)[1:-1],
                    np.maximum(np.asarray(f_zt)[1:-1],
                               np.asarray(f_zt)[2:]))
    assert (tm[1:-1] >= lo - 1e-12).all()
    assert (tm[1:-1] <= hi + 1e-12).all()
    assert tm[0] == np.asarray(f_zm)[0]
    assert tm[-1] == np.asarray(f_zm)[-1]


# --------------------------------------- Tier-0: clip_rcm / cloud cover
def test_tier0_clip_rcm():
    rtm = jnp.asarray([1e-3, 1e-3, 5e-6, 0.0, 2e-2])
    rcm = jnp.asarray([2e-4, 1e-3, 1e-5, 1e-6, 2.5e-2])
    out = np.asarray(pcd.clip_rcm(rtm, rcm))
    # untouched where rcm <= rtm (bitwise)
    assert out[0] == 2e-4 and out[1] == 1e-3
    # clipped to rtm - epsilon (never negative)
    np.testing.assert_allclose(out[2], 5e-6 - EPSILON_R8, rtol=0)
    assert out[3] == 0.0
    np.testing.assert_allclose(out[4], 2e-2 - EPSILON_R8, rtol=0)
    assert (out >= 0.0).all()


def test_tier0_cloud_cover_properties(grid_f):
    gr = grid_f
    nz = gr.zt.shape[0]
    rng = np.random.default_rng(4)
    mixt_frac = jnp.asarray(rng.uniform(0.2, 0.8, nz))
    chi_1 = jnp.asarray(rng.normal(0.0, 1e-4, nz))
    chi_2 = jnp.asarray(rng.normal(0.0, 1e-4, nz))

    # clear column: everything zero in -> everything zero out
    zero = jnp.zeros(nz)
    cc, ril = pcd.compute_cloud_cover(gr, mixt_frac, chi_1, chi_2,
                                      zero, zero)
    assert np.asarray(cc).max() == 0.0 and np.asarray(ril).max() == 0.0

    # a broken cloud profile: cloud cover boosts cloud_frac, in-layer
    # rcm boosts rcm, and everything stays realizable
    rcm = np.zeros(nz)
    rcm[10:14] = 3e-5
    rcm[20] = 2e-5   # one-level cloud (both top and base branches)
    cf = np.clip(rcm * 2e4, 0.0, 1.0)
    cc, ril = pcd.compute_cloud_cover(
        gr, mixt_frac, chi_1, chi_2, jnp.asarray(cf), jnp.asarray(rcm))
    cc, ril = np.asarray(cc), np.asarray(ril)
    assert (cc >= cf - 1e-15).all()          # cc = cf / vcf, vcf <= 1
    assert (ril >= rcm - 1e-20).all()
    assert (cc <= 1.0 + 1e-12).all() and (ril >= 0.0).all()
    # away from cloud boundaries nothing changes
    interior = (rcm >= RC_TOL)
    interior[1:] &= (rcm[:-1] >= RC_TOL)
    interior[:-1] &= (rcm[1:] >= RC_TOL)
    np.testing.assert_array_equal(cc[interior], cf[interior])
    # boundary levels of the partial cloud are boosted
    assert cc[10] > cf[10] and cc[13] > cf[13] and cc[20] > cf[20]


# ------------------------------------ Tier-0: driver-level realizability
def test_tier0_driver_realizability(gdrv, tunables, replay_syn,
                                    replay_adv):
    for outs, _, _ in (replay_syn, replay_adv):
        for name in ["cloud_frac", "cloud_cover", "cloud_frac_zm",
                     "ice_supersat_frac", "ice_supersat_frac_zm"]:
            v = outs[name]
            assert (v >= 0.0).all() and (v <= 1.0).all(), name
        for name in ["rcm", "rcm_in_layer", "rcm_zm", "rcp2",
                     "rcp2_zt", "sigma_sqd_w", "sigma_sqd_w_zt"]:
            assert (outs[name] >= 0.0).all(), name
        # l_use_cloud_cover substitution: rcm == rcm_in_layer and
        # cloud_frac == min(1, cloud_cover)
        np.testing.assert_array_equal(outs["rcm"],
                                      outs["rcm_in_layer"])
        np.testing.assert_array_equal(
            outs["cloud_frac"], np.minimum(1.0, outs["cloud_cover"]))
        assert (outs["rcm_supersat_adj"] == 0.0).all()


def test_tier0_driver_clear_column(grid_f, tunables):
    """A bone-dry column must produce zero cloud everywhere through
    the full driver (both pdf calls, trapezoids and cloud cover)."""
    tun = tunables[0]
    gr = grid_f
    nz = gr.zt.shape[0]
    z_t = np.maximum(np.asarray(gr.zt), 0.0)
    p = 1.0e5 * np.exp(-z_t / 8000.0)
    p[0] = p[1]
    exner = (p / 1.0e5) ** (287.042 / 1004.64)
    T = np.maximum(290.0 - 0.006 * z_t, 200.0)
    thlm = T / exner
    from eam_jax.clubb import saturation as csat
    rsl = np.asarray(csat.sat_mixrat_liq(p, T))
    rtm = np.clip(0.2 * rsl, 1e-7, 0.02)   # 20% RH
    zer = np.zeros(nz)
    wp2 = np.full(nz, 0.3)
    wp3 = np.full(nz, 0.05)
    rtp2 = (0.05 * rtm) ** 2
    thlp2 = np.full(nz, 0.09)
    args = [0.3 * np.sqrt(wp2 * rtp2), thlm, -0.2 * np.sqrt(wp2 * thlp2),
            rtp2, 0.1 * rtp2 ** 1.5, thlp2, zer, np.zeros(nz), wp2, wp3,
            zer, zer, np.full(nz, 5.0), np.full(nz, 0.2),
            np.full(nz, -0.02), np.full(nz, -2.0), np.full(nz, 0.2),
            np.full(nz, 0.01), p, exner, thlm, thlm, rtm]
    outs, pz, pm = pcd.pdf_closure_driver(
        gr, *[jnp.asarray(a) for a in args], **tun)
    for name in ["rcm", "cloud_frac", "cloud_cover", "rcm_in_layer",
                 "cloud_frac_zm", "rcm_zm", "wprcp", "rcp2", "rcp2_zt"]:
        assert np.abs(np.asarray(outs[name])).max() == 0.0, name
