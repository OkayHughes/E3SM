"""Tier-1 golden replay + Tier-0 properties for eam_jax.clubb
(slice 1: grid, saturation, tridiag/dgtsv, pdf_closure ADG1)."""

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
from eam_jax.clubb import saturation as csat  # noqa: E402
from eam_jax.clubb import tridiag as ctri  # noqa: E402

GOLDEN = ROOT / "golden"

MOMENT_SLOTS = ["wp2rtp", "wp2thlp", "cloud_frac", "ice_supersat_frac",
                "rcm", "wpthvp", "wp2thvp", "rtpthvp", "thlpthvp",
                "wprcp", "wp2rcp", "rtprcp", "thlprcp", "rcp2",
                "uprcp", "vprcp", "rc_coef"]

PDFP_SLOTS = ["w_1", "w_2", "varnce_w_1", "varnce_w_2", "rt_1", "rt_2",
              "varnce_rt_1", "varnce_rt_2", "thl_1", "thl_2",
              "varnce_thl_1", "varnce_thl_2", "corr_w_rt_1",
              "corr_w_rt_2", "corr_w_thl_1", "corr_w_thl_2",
              "corr_rt_thl_1", "corr_rt_thl_2", "alpha_thl", "alpha_rt",
              "crt_1", "crt_2", "cthl_1", "cthl_2", "chi_1", "chi_2",
              "stdev_chi_1", "stdev_chi_2", "stdev_eta_1", "stdev_eta_2",
              "covar_chi_eta_1", "covar_chi_eta_2", "corr_w_chi_1",
              "corr_w_chi_2", "corr_w_eta_1", "corr_w_eta_2",
              "corr_chi_eta_1", "corr_chi_eta_2", "rsatl_1", "rsatl_2",
              "rc_1", "rc_2", "cloud_frac_1", "cloud_frac_2",
              "mixt_frac", "ice_supersat_frac_1", "ice_supersat_frac_2"]

# MEASURED agreement (40 cols x 73 levs): every output that does not
# pass through erf()/exp() replays BITWISE (all PDF-component means,
# variances, correlations, chi/eta transform, wp2rtp/wp2thlp, rc_coef,
# mixt_frac).  Only the partly-cloudy-branch quantities differ, by
# cross-libm (gfortran vs XLA) erf/exp tails: cloud_frac abs <= 2.3e-16
# (1 ulp of 1), rc_i abs <= 1.7e-19 (rel 4.5e-9 at |rc| ~ 2e-11
# leading edges), rcp2 abs <= 1.8e-22, and their downstream x'rc'/th_v
# moments at the same absolute level (wp2thvp/wp2rcp "rel" errors are
# sign flips of ~1e-300 cancellation residue).  The atol floors below
# are ~10x the measured maxima; rtol stays at the 1e-12 target.
MOMENT_ATOL = {"wp2rtp": 0.0, "wp2thlp": 0.0, "rc_coef": 0.0,
               "cloud_frac": 2e-15, "ice_supersat_frac": 2e-15,
               "rcm": 1e-18, "wpthvp": 2e-15, "wp2thvp": 2e-15,
               "rtpthvp": 2e-17, "thlpthvp": 2e-15, "wprcp": 1e-18,
               "wp2rcp": 1e-18, "rtprcp": 1e-20, "thlprcp": 1e-18,
               "rcp2": 1e-20, "uprcp": 1e-18, "vprcp": 1e-18}
PDFP_ATOL = {"cloud_frac_1": 2e-15, "cloud_frac_2": 2e-15,
             "rc_1": 1e-18, "rc_2": 1e-18,
             "ice_supersat_frac_1": 2e-15, "ice_supersat_frac_2": 2e-15}


def _load(name):
    path = GOLDEN / name
    if not path.exists():
        pytest.skip("golden archive missing (run harness in container)")
    return np.load(path)


@pytest.fixture(scope="module")
def ggrid():
    return _load("clubb_grid.npz")


@pytest.fixture(scope="module")
def gtri():
    return _load("clubb_tridag.npz")


@pytest.fixture(scope="module")
def gsat():
    return _load("clubb_sat.npz")


@pytest.fixture(scope="module")
def gpdf():
    return _load("clubb_pdf_closure.npz")


# ---------------------------------------------------------------- grid
@pytest.mark.parametrize("tag", ["a", "b", "c"])
def test_grid_arrays_golden(ggrid, tag):
    gr = cgrid.setup_grid(ggrid[f"{tag}_zi_in"], ggrid[f"{tag}_zt_in"])
    for name, ours in [("zm", gr.zm), ("zt", gr.zt), ("dzm", gr.dzm),
                       ("dzt", gr.dzt), ("invrs_dzm", gr.invrs_dzm),
                       ("invrs_dzt", gr.invrs_dzt)]:
        np.testing.assert_allclose(np.asarray(ours),
                                   ggrid[f"{tag}_{name}"], rtol=1e-14,
                                   err_msg=name)
    w_zt2zm = ggrid[f"{tag}_w_zt2zm"]
    w_zm2zt = ggrid[f"{tag}_w_zm2zt"]
    np.testing.assert_allclose(np.asarray(gr.factor_zt2zm), w_zt2zm[0],
                               rtol=1e-14)
    np.testing.assert_allclose(1.0 - np.asarray(gr.factor_zt2zm),
                               w_zt2zm[1], rtol=1e-13, atol=1e-15)
    np.testing.assert_allclose(np.asarray(gr.factor_zm2zt), w_zm2zt[0],
                               rtol=1e-14)
    np.testing.assert_allclose(1.0 - np.asarray(gr.factor_zm2zt),
                               w_zm2zt[1], rtol=1e-13, atol=1e-15)


@pytest.mark.parametrize("tag", ["a", "b", "c"])
def test_grid_ops_golden(ggrid, tag):
    gr = cgrid.setup_grid(ggrid[f"{tag}_zi_in"], ggrid[f"{tag}_zt_in"])
    fzt = ggrid[f"{tag}_fields_zt"]
    fzm = ggrid[f"{tag}_fields_zm"]
    for i in range(fzt.shape[0]):
        np.testing.assert_allclose(np.asarray(cgrid.zt2zm(gr, fzt[i])),
                                   ggrid[f"{tag}_zt2zm"][i], rtol=1e-13,
                                   atol=1e-15)
        np.testing.assert_allclose(np.asarray(cgrid.zm2zt(gr, fzm[i])),
                                   ggrid[f"{tag}_zm2zt"][i], rtol=1e-13,
                                   atol=1e-15)
        np.testing.assert_allclose(np.asarray(cgrid.ddzt(gr, fzt[i])),
                                   ggrid[f"{tag}_ddzt"][i], rtol=1e-13,
                                   atol=1e-16)
        np.testing.assert_allclose(np.asarray(cgrid.ddzm(gr, fzm[i])),
                                   ggrid[f"{tag}_ddzm"][i], rtol=1e-13,
                                   atol=1e-16)


def test_grid_ops_properties(ggrid):
    gr = cgrid.setup_grid(ggrid["a_zi_in"], ggrid["a_zt_in"])
    const = jnp.full(gr.zt.shape, 3.5)
    np.testing.assert_allclose(np.asarray(cgrid.zt2zm(gr, const)), 3.5,
                               rtol=1e-15)
    np.testing.assert_allclose(np.asarray(cgrid.zm2zt(gr, const)), 3.5,
                               rtol=1e-15)
    np.testing.assert_allclose(np.asarray(cgrid.ddzt(gr, const)), 0.0,
                               atol=1e-18)
    np.testing.assert_allclose(np.asarray(cgrid.ddzm(gr, const)), 0.0,
                               atol=1e-18)
    # derivative of a linear-in-z field is its slope
    lin_t = 2.0 * gr.zt + 1.0
    np.testing.assert_allclose(np.asarray(cgrid.ddzt(gr, lin_t)), 2.0,
                               rtol=1e-12)
    lin_m = -0.5 * gr.zm + 4.0
    np.testing.assert_allclose(np.asarray(cgrid.ddzm(gr, lin_m)), -0.5,
                               rtol=1e-12)


# -------------------------------------------------------------- tridag
def test_tridag_golden(gtri):
    """Measured envelope (see PORT_NOTES in tridiag.py): most cases
    replay bitwise (12/18, all 73-level CLUBB-like dominant systems),
    the rest differ by single ulps because Ubuntu's prebuilt liblapack
    was compiled with FMA contraction (aarch64 gfortran default) —
    max 4.7e-14 elementwise on well-conditioned systems.  The
    deliberately ill-conditioned mixed-magnitude family
    (case index %3 == 2) amplifies that single-ulp noise to 5.2e-11
    solution-scaled; there the assertion is on the scaled error."""
    ncases = int(gtri["ncases"])
    for i in range(ncases):
        sol, sing = ctri.tridag_solve(gtri[f"{i}_supd"], gtri[f"{i}_diag"],
                                      gtri[f"{i}_subd"], gtri[f"{i}_rhs"])
        assert not bool(sing)
        ref = gtri[f"{i}_solution"]
        sol = np.asarray(sol)
        scale = np.abs(ref).max()
        if i % 3 == 2:  # ill-conditioned family
            assert np.abs(sol - ref).max() / scale < 1e-9, f"case {i}"
        else:
            np.testing.assert_allclose(sol, ref, rtol=1e-12,
                                       atol=1e-15 * scale,
                                       err_msg=f"case {i}")


def test_tridag_residual_and_singular():
    rng = np.random.default_rng(3)
    n = 30
    sub = rng.normal(size=n)
    dia = rng.normal(size=n) * 0.2
    sup = rng.normal(size=n)
    rhs = rng.normal(size=(n, 2))
    sol, sing = ctri.tridag_solve(sup, dia, sub, rhs)
    A = np.diag(dia) + np.diag(sub[1:], -1) + np.diag(sup[:-1], 1)
    assert not bool(sing)
    np.testing.assert_allclose(A @ np.asarray(sol), rhs, atol=1e-10)
    # singular system flagged, CLUBB's -999 fill applied
    z = np.zeros(4)
    sol, sing = ctri.tridag_solve(z, z, z, np.ones(4))
    assert bool(sing)
    assert np.all(np.asarray(sol) == -999.0)


# ----------------------------------------------------------------- sat
def test_sat_golden(gsat):
    rsl = np.asarray(csat.sat_mixrat_liq(gsat["p"], gsat["t"]))
    rsi = np.asarray(csat.sat_mixrat_ice(gsat["p"], gsat["t"]))
    np.testing.assert_allclose(rsl, gsat["rsl"], rtol=1e-13)
    np.testing.assert_allclose(rsi, gsat["rsi"], rtol=1e-13)


def test_sat_properties():
    t = np.linspace(180.0, 305.0, 200)
    p = np.full_like(t, 8.0e4)
    rsl = np.asarray(csat.sat_mixrat_liq(p, t))
    rsi = np.asarray(csat.sat_mixrat_ice(p, t))
    assert (rsl > 0).all() and (rsi > 0).all()
    # monotonically increasing in T away from the polynomial clip
    assert (np.diff(rsl[t > 190.0]) > 0).all()
    # ice SVP below liquid SVP under freezing (within fit validity)
    m = (t > 200.0) & (t < 273.0)
    assert (rsi[m] < rsl[m]).all()


# --------------------------------------------------------- pdf_closure
def _pdf_batch(gpdf):
    meta = json.loads(str(gpdf["meta"]))
    params = gpdf["params"]
    idx = meta["param_indices"]
    beta = float(params[idx["beta"]])
    mf_max = pc.mixt_frac_max_mag(float(params[idx["Skw_max_mag"]]))
    ins = [gpdf[f"in_{k}"] for k in
           ["p", "exner", "thv_ds", "wm", "wp2", "wp3", "sigma_sqd_w",
            "skw", "skthl", "skrt", "rtm", "rtp2", "wprtp", "thlm",
            "thlp2", "wpthlp", "um", "up2", "upwp", "vm", "vp2", "vpwp",
            "rtpthlp"]]
    fn = jax.vmap(lambda *a: pc.pdf_closure(*a, beta=beta,
                                            mf_max_mag=mf_max))
    return fn(*[jnp.asarray(x) for x in ins])


@pytest.fixture(scope="module")
def pdf_out(gpdf):
    moments, pdfp = _pdf_batch(gpdf)
    return ({k: np.asarray(v) for k, v in moments.items()},
            {k: np.asarray(v) for k, v in pdfp.items()})


def test_pdf_closure_moments_golden(gpdf, pdf_out):
    moments, _ = pdf_out
    ref = gpdf["moments"]
    for j, name in enumerate(MOMENT_SLOTS):
        np.testing.assert_allclose(
            moments[name], ref[:, :, j], rtol=1e-12,
            atol=MOMENT_ATOL.get(name, 0.0), err_msg=name)


def test_pdf_closure_pdf_params_golden(gpdf, pdf_out):
    _, pdfp = pdf_out
    ref = gpdf["pdf_params"]
    for j, name in enumerate(PDFP_SLOTS):
        np.testing.assert_allclose(
            pdfp[name], ref[:, :, j], rtol=1e-12,
            atol=PDFP_ATOL.get(name, 0.0), err_msg=name)


def test_pdf_closure_sigma_unchanged(gpdf):
    # ADG1 does not overwrite sigma_sqd_w (only ADG2 does)
    np.testing.assert_array_equal(gpdf["sigma_sqd_w_out"],
                                  gpdf["in_sigma_sqd_w"])


# ------------------------------------------------- Tier-0: realizability
def test_tier0_realizability(gpdf, pdf_out):
    meta = json.loads(str(gpdf["meta"]))
    params = gpdf["params"]
    mf_max = pc.mixt_frac_max_mag(
        float(params[meta["param_indices"]["Skw_max_mag"]]))
    moments, pdfp = pdf_out
    for v in ["varnce_w_1", "varnce_w_2", "varnce_rt_1", "varnce_rt_2",
              "varnce_thl_1", "varnce_thl_2"]:
        assert (pdfp[v] >= 0.0).all(), v
    for v in ["corr_rt_thl_1", "corr_rt_thl_2", "corr_chi_eta_1",
              "corr_chi_eta_2"]:
        assert (np.abs(pdfp[v]) <= 0.99).all(), v
    for v in ["cloud_frac_1", "cloud_frac_2", "ice_supersat_frac_1",
              "ice_supersat_frac_2"]:
        assert ((pdfp[v] >= 0.0) & (pdfp[v] <= 1.0)).all(), v
    assert ((moments["cloud_frac"] >= 0.0)
            & (moments["cloud_frac"] <= 1.0)).all()
    assert ((moments["ice_supersat_frac"] >= 0.0)
            & (moments["ice_supersat_frac"] <= 1.0)).all()
    assert (moments["rcm"] >= 0.0).all()
    assert (moments["rcp2"] >= 0.0).all()
    assert (pdfp["rc_1"] >= 0.0).all() and (pdfp["rc_2"] >= 0.0).all()
    mf = pdfp["mixt_frac"]
    assert ((mf >= 1.0 - mf_max) & (mf <= mf_max)).all()
    assert (pdfp["stdev_chi_1"] >= 0.0).all()
    assert (pdfp["stdev_eta_1"] >= 0.0).all()


def _column_case(rh, skw_val, wprtp_corr=0.7, t_sfc=295.0, nz=40):
    """Build a simple consistent single column for Tier-0 cases."""
    z = np.linspace(50.0, 12000.0, nz)
    p = 1.0e5 * np.exp(-z / 8000.0)
    exner = (p / 1.0e5) ** (287.042 / 1004.64)
    T = np.maximum(t_sfc - 0.0065 * z, 200.0)
    thlm = T / exner
    rsl = np.asarray(csat.sat_mixrat_liq(p, T))
    rtm = np.clip(rh * rsl, 1e-7, 0.03)
    wp2 = np.full(nz, 0.36)
    wp3 = skw_val * wp2 ** 1.5
    rtp2 = (0.1 * rtm) ** 2
    thlp2 = np.full(nz, 0.25)
    wprtp = wprtp_corr * np.sqrt(wp2 * rtp2)
    wpthlp = 0.3 * np.sqrt(wp2 * thlp2)
    rtpthlp = -0.4 * np.sqrt(rtp2 * thlp2)
    up2 = np.full(nz, 0.2)
    vp2 = np.full(nz, 0.2)
    upwp = np.full(nz, -0.02)
    vpwp = np.full(nz, 0.01)
    skw = np.asarray(pc.skx_func(jnp.asarray(wp2), jnp.asarray(wp3),
                                 2.0e-2, 0.0))
    gam = np.asarray(pc.gamma_skw_fnc(jnp.asarray(skw), 0.12, 0.28, 1.2))
    sig = np.asarray(pc.compute_sigma_sqd_w(
        jnp.asarray(gam), jnp.asarray(wp2), jnp.asarray(thlp2),
        jnp.asarray(rtp2), jnp.asarray(up2), jnp.asarray(vp2),
        jnp.asarray(wpthlp), jnp.asarray(wprtp), jnp.asarray(upwp),
        jnp.asarray(vpwp)))
    args = [p, exner, thlm.copy(), np.zeros(nz), wp2, wp3, sig, skw,
            np.zeros(nz), np.zeros(nz), rtm, rtp2, wprtp, thlm, thlp2,
            wpthlp, np.full(nz, 5.0), up2, upwp, np.full(nz, -2.0),
            vp2, vpwp, rtpthlp]
    beta = 2.4
    mf_max = pc.mixt_frac_max_mag(4.5)
    return pc.pdf_closure(*[jnp.asarray(a) for a in args],
                          beta=beta, mf_max_mag=mf_max)


def test_tier0_clear_and_saturated_limits():
    m_dry, p_dry = _column_case(rh=0.25, skw_val=0.8)
    assert np.asarray(m_dry["cloud_frac"]).max() == 0.0
    assert np.asarray(m_dry["rcm"]).max() == 0.0
    # no cloud => zero liquid-water fluxes and rcp2
    assert np.abs(np.asarray(m_dry["wprcp"])).max() == 0.0
    assert np.asarray(m_dry["rcp2"]).max() == 0.0

    m_sat, p_sat = _column_case(rh=1.5, skw_val=0.5)
    cf = np.asarray(m_sat["cloud_frac"])
    warm = np.asarray(p_sat["rsatl_1"]) > 1e-3  # lower troposphere
    assert cf[warm].min() > 0.99
    assert np.asarray(m_sat["rcm"])[warm].min() > 0.0


def test_tier0_convective_flux_signs():
    # moist updrafts (wprtp > 0) in a partly-cloudy convective column
    # must carry positive liquid-water flux and enhance buoyancy flux
    m, p = _column_case(rh=0.97, skw_val=1.5, wprtp_corr=0.8)
    cloudy = np.asarray(m["cloud_frac"]) > 0.05
    assert cloudy.any()
    assert (np.asarray(m["wprcp"])[cloudy] > 0.0).all()
    wpthvp = np.asarray(m["wpthvp"])
    wpthlp_plus_vap = np.asarray(m["wpthvp"] - m["rc_coef"] * m["wprcp"])
    assert (wpthvp[cloudy] > wpthlp_plus_vap[cloudy]).all()


def test_tier0_mixt_frac_skewness_symmetry():
    # wm = 0: negating Skw mirrors the PDF: mixt_frac -> 1 - mixt_frac,
    # w_1 -> -w_2 (components swap)
    m_p, p_p = _column_case(rh=0.5, skw_val=2.0)
    m_n, p_n = _column_case(rh=0.5, skw_val=-2.0)
    np.testing.assert_allclose(np.asarray(p_n["mixt_frac"]),
                               1.0 - np.asarray(p_p["mixt_frac"]),
                               rtol=1e-12)
    np.testing.assert_allclose(np.asarray(p_n["w_1"]),
                               -np.asarray(p_p["w_2"]), rtol=1e-12)


def test_tier0_zero_variance_w_degenerate():
    # wp2 below w_tol_sqd: single-Gaussian responders, w components
    # collapse onto the mean
    nz = 10
    p = np.full(nz, 8.0e4)
    exner = (p / 1.0e5) ** (287.042 / 1004.64)
    T = np.full(nz, 285.0)
    thlm = T / exner
    rtm = np.full(nz, 0.008)
    zero = np.zeros(nz)
    args = [p, exner, thlm, zero, np.full(nz, 1e-6), zero,
            np.full(nz, 0.2), zero, zero, zero, rtm,
            np.full(nz, 1e-8), zero, thlm, np.full(nz, 0.04), zero,
            np.full(nz, 3.0), np.full(nz, 0.1), zero,
            np.full(nz, 1.0), np.full(nz, 0.1), zero, zero]
    m, pp = pc.pdf_closure(*[jnp.asarray(a) for a in args], beta=2.4,
                           mf_max_mag=pc.mixt_frac_max_mag(4.5))
    np.testing.assert_array_equal(np.asarray(pp["mixt_frac"]), 0.5)
    np.testing.assert_array_equal(np.asarray(pp["w_1"]), 0.0)
    np.testing.assert_array_equal(np.asarray(pp["varnce_w_1"]), 0.0)
    np.testing.assert_array_equal(np.asarray(pp["rt_1"]),
                                  np.asarray(pp["rt_2"]))
    # single-Gaussian responder keeps the full variance
    np.testing.assert_array_equal(np.asarray(pp["varnce_rt_1"]),
                                  np.full(nz, 1e-8))
