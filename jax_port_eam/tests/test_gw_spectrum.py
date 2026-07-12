"""Tier-1 golden replay + Tier-0 properties for the full-spectrum
(non-orographic) gravity-wave suite (eam_jax.gw_spectrum)."""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
GOLDEN = ROOT / "golden" / "gw_spectrum_golden.npz"

from eam_jax import gw_spectrum as gws  # noqa: E402

RTOL = 1e-12


@pytest.fixture(scope="module")
def gold():
    if not GOLDEN.exists():
        pytest.skip("golden archive missing (run harness in container)")
    return np.load(GOLDEN)


@pytest.fixture(scope="module")
def meta(gold):
    return json.loads(str(gold["__metadata__"]))


@pytest.fixture(scope="module")
def params(gold, meta):
    p = meta["params"]
    return gws.make_gw_spectrum_params(
        gold["alpha"], p["pgwv"], p["dc"], p["kbotbg"], meta["nlev"],
        fcrit2=p["fcrit2"], kwv=p["kwv"], gravit=p["gravit"],
        rair=p["rair"], ktop=p["ktop"], tau_0_ubc=p["tau_0_ubc"])


@pytest.fixture(scope="module")
def fp(gold, meta, params):
    p = meta["params"]
    return gws.make_front_params(params, p["taubgnd"], p["frontgfc"],
                                 p["kfront"])


def beres_params(gold, meta, params, use_old):
    p = meta["params"]
    return gws.make_beres_params(
        params, gold["mfcc"], gold["pref_edge"],
        plev_src_wind=p["gw_convect_plev_src_wind"],
        maxq0_conversion_factor=p["gw_convect_hcf"],
        hdepth_scaling_factor=p["hdepth_scaling_factor"],
        hdepth_min=p["gw_convect_hdepth_min"],
        storm_speed_min=p["gw_convect_storm_speed_min"],
        use_gw_convect_old=use_old)


def run_beres(gold, meta, params, use_old):
    bp = beres_params(gold, meta, params, use_old)
    return gws.gw_beres_src(gold["lat"], gold["u"], gold["v"],
                            gold["netdt"], gold["zm"], params, bp)


def run_drag(gold, meta, params, src, effgw, do_taper):
    """src = (src_level, tend_level, tau, ubm, ubi, xv, yv, c)."""
    src_level, tend_level, tau, ubm, ubi, xv, yv, c = src
    return gws.gw_drag_prof(
        src_level, tend_level, meta["dt"], gold["lat"], gold["t"],
        gold["ti"], gold["pmid"], gold["pint"], gold["dpm"],
        gold["rdpm"], gold["piln"], gold["rhoi"], gold["nm"],
        gold["ni"], ubm, ubi, xv, yv, effgw, c, gold["q"], gold["dse"],
        tau, params, do_taper=do_taper)


# ---------------------------------------------------------------------------
# Tier-1 golden replay
# ---------------------------------------------------------------------------
def test_k_src_wind_matches_init(gold, meta, params):
    bp = beres_params(gold, meta, params, True)
    assert bp["k_src_wind"] == meta["params"]["k_src_wind"]


@pytest.mark.parametrize("use_old,sfx", [(True, "_old"), (False, "_new")])
def test_gw_beres_src(gold, meta, params, use_old, sfx):
    (src, tend, tau, ubm, ubi, xv, yv, c, hdepth,
     maxq0) = run_beres(gold, meta, params, use_old)
    np.testing.assert_array_equal(np.asarray(src),
                                  gold["bsrc_level" + sfx])
    np.testing.assert_array_equal(np.asarray(tend),
                                  gold["btend_level" + sfx])
    np.testing.assert_allclose(np.asarray(tau), gold["btau" + sfx],
                               rtol=RTOL, atol=1e-300)
    np.testing.assert_allclose(np.asarray(ubm), gold["bubm" + sfx],
                               rtol=RTOL, atol=1e-14)
    np.testing.assert_allclose(np.asarray(ubi), gold["bubi" + sfx],
                               rtol=RTOL, atol=1e-14)
    np.testing.assert_allclose(np.asarray(xv), gold["bxv" + sfx],
                               rtol=1e-13)
    np.testing.assert_allclose(np.asarray(yv), gold["byv" + sfx],
                               rtol=1e-13)
    np.testing.assert_array_equal(np.asarray(c), gold["bc" + sfx])
    np.testing.assert_allclose(np.asarray(hdepth),
                               gold["bhdepth" + sfx], rtol=RTOL)
    np.testing.assert_allclose(np.asarray(maxq0), gold["bmaxq0" + sfx],
                               rtol=RTOL, atol=1e-300)


def test_gw_cm_src(gold, meta, params, fp):
    src, tend, tau, ubm, ubi, xv, yv, c = gws.gw_cm_src(
        gold["u"], gold["v"], gold["frontgf"], params, fp)
    np.testing.assert_array_equal(np.asarray(src), gold["csrc_level"])
    np.testing.assert_array_equal(np.asarray(tend),
                                  gold["ctend_level"])
    np.testing.assert_allclose(np.asarray(tau), gold["ctau"],
                               rtol=RTOL, atol=1e-300)
    np.testing.assert_allclose(np.asarray(ubm), gold["cubm"],
                               rtol=RTOL, atol=1e-14)
    np.testing.assert_allclose(np.asarray(ubi), gold["cubi"],
                               rtol=RTOL, atol=1e-14)
    np.testing.assert_allclose(np.asarray(xv), gold["cxv"], rtol=1e-13)
    np.testing.assert_allclose(np.asarray(yv), gold["cyv"], rtol=1e-13)
    np.testing.assert_allclose(np.asarray(c), gold["cc"], rtol=RTOL)


def test_gwd_project_tau(gold, params):
    taucd = gws._project_tau(gold["btend_level_old"],
                             gold["btau_old"], gold["bubi_old"],
                             gold["bc_old"], gold["bxv_old"],
                             gold["byv_old"], params)
    np.testing.assert_allclose(np.asarray(taucd), gold["ptaucd"],
                               rtol=RTOL, atol=1e-300)


DRAG_RUNS = [("bd", "b", "_old", "effgw_beres", False),
             ("cd", "c", "", "effgw_cm", False),
             ("cdt", "c", "", "effgw_cm", True)]


@pytest.mark.parametrize("tag,src,sfx,effkey,taper", DRAG_RUNS)
def test_gw_drag_prof(gold, meta, params, tag, src, sfx, effkey, taper):
    inp = tuple(gold[f"{src}{k}{sfx}"] for k in
                ("src_level", "tend_level", "tau", "ubm", "ubi", "xv",
                 "yv", "c"))
    out = run_drag(gold, meta, params, inp, meta["params"][effkey],
                   taper)
    names = ("tau", "utgw", "vtgw", "ttgw", "qtgw", "taucd", "egwdffi",
             "gwut", "dttdf", "dttke")
    atols = dict(tau=1e-300, utgw=1e-20, vtgw=1e-20, ttgw=1e-18,
                 qtgw=1e-25, taucd=1e-300, egwdffi=1e-300, gwut=1e-20,
                 dttdf=1e-18, dttke=1e-18)
    for name, got in zip(names, out):
        np.testing.assert_allclose(
            np.asarray(got), gold[f"{tag}_{name}"], rtol=RTOL,
            atol=atols[name], err_msg=f"{tag}_{name}")


def test_vd_lu(gold, meta):
    p = meta["params"]
    ntop, nbot = p["lu_ntop"], p["lu_nbot"]
    decomp = gws.vd_lu_decomp(gold["lu_ksrf"], gold["lu_kv"],
                              gold["lu_tmpi"], gold["rdpm"],
                              meta["dt"], p["gravit"],
                              gold["lu_cc_top"], ntop, nbot)
    for k, name in (("ca", "lu_ca"), ("cc", "lu_cc"),
                    ("dnom", "lu_dnom"), ("ze", "lu_ze")):
        np.testing.assert_allclose(np.asarray(decomp[k]), gold[name],
                                   rtol=RTOL, atol=1e-300, err_msg=k)
    import jax.numpy as jnp
    for m in range(gold["lu_q"].shape[2]):
        qs = gws.vd_lu_solve(jnp.asarray(gold["lu_q"][:, :, m]),
                             decomp, ntop, nbot, gold["lu_cd_top"])
        np.testing.assert_allclose(np.asarray(qs),
                                   gold["lu_q_out"][:, :, m],
                                   rtol=RTOL, atol=1e-300,
                                   err_msg=f"lu q{m}")


def test_momentum_energy_conservation(gold, meta, params):
    import jax.numpy as jnp
    r = gws.momentum_energy_conservation(
        gold["btend_level_old"], meta["dt"], gold["bd_taucd"],
        gold["pint"], gold["dpm"], gold["u"], gold["v"],
        jnp.asarray(gold["bd_utgw"]), jnp.asarray(gold["bd_vtgw"]),
        jnp.asarray(gold["bd_ttgw"]), jnp.asarray(gold["bd_utgw"]),
        jnp.asarray(gold["bd_vtgw"]), jnp.asarray(gold["bd_ttgw"]),
        params)
    names = ("mec_b_dudt", "mec_b_dvdt", "mec_b_dsdt", "mec_b_utgw",
             "mec_b_vtgw", "mec_b_ttgw")
    for name, got in zip(names, r):
        np.testing.assert_allclose(np.asarray(got), gold[name],
                                   rtol=RTOL, atol=1e-18,
                                   err_msg=name)


def test_chained_replay(gold, meta, params, fp):
    """gw_tend order: prof -> beres src -> drag -> conservation ->
    cm src -> drag -> conservation, accumulating ptend."""
    import jax.numpy as jnp
    p = meta["params"]
    rhoi, ti, nm, ni = gws.gw_prof(p["cpair"], gold["t"], gold["pmid"],
                                   gold["pint"], params)
    np.testing.assert_allclose(np.asarray(rhoi), gold["rhoi"],
                               rtol=1e-13)
    np.testing.assert_allclose(np.asarray(ni), gold["ni"], rtol=1e-13)

    bsrc = run_beres(gold, meta, params, True)
    (src_b, tend_b, tau_b, ubm_b, ubi_b, xv_b, yv_b, c_b, _, _) = bsrc
    out_b = gws.gw_drag_prof(
        src_b, tend_b, meta["dt"], gold["lat"], gold["t"], ti,
        gold["pmid"], gold["pint"], gold["dpm"], gold["rdpm"],
        gold["piln"], rhoi, nm, ni, ubm_b, ubi_b, xv_b, yv_b,
        p["effgw_beres"], c_b, gold["q"], gold["dse"], tau_b, params)
    (_, utgw, vtgw, ttgw, qtgw, taucd, _, _, _, _) = out_b
    ptu, ptv, pts, utgw, vtgw, ttgw = gws.momentum_energy_conservation(
        tend_b, meta["dt"], taucd, gold["pint"], gold["dpm"],
        gold["u"], gold["v"], utgw, vtgw, ttgw,
        jnp.array(utgw), jnp.array(vtgw), jnp.array(ttgw), params)
    ptq = qtgw

    csrc = gws.gw_cm_src(gold["u"], gold["v"], gold["frontgf"], params,
                         fp)
    (src_c, tend_c, tau_c, ubm_c, ubi_c, xv_c, yv_c, c_c) = csrc
    out_c = gws.gw_drag_prof(
        src_c, tend_c, meta["dt"], gold["lat"], gold["t"], ti,
        gold["pmid"], gold["pint"], gold["dpm"], gold["rdpm"],
        gold["piln"], rhoi, nm, ni, ubm_c, ubi_c, xv_c, yv_c,
        p["effgw_cm"], c_c, gold["q"], gold["dse"], tau_c, params)
    (_, utgw_c, vtgw_c, ttgw_c, qtgw_c, taucd_c, _, _, _, _) = out_c
    ptu = ptu + utgw_c
    ptv = ptv + vtgw_c
    pts = pts + ttgw_c
    ptq = ptq + qtgw_c
    ptu, ptv, pts, _, _, _ = gws.momentum_energy_conservation(
        tend_c, meta["dt"], taucd_c, gold["pint"], gold["dpm"],
        gold["u"], gold["v"], ptu, ptv, pts, utgw_c, vtgw_c, ttgw_c,
        params)

    np.testing.assert_allclose(np.asarray(ptu), gold["chain_dudt"],
                               rtol=RTOL, atol=1e-18)
    np.testing.assert_allclose(np.asarray(ptv), gold["chain_dvdt"],
                               rtol=RTOL, atol=1e-18)
    np.testing.assert_allclose(np.asarray(pts), gold["chain_dsdt"],
                               rtol=RTOL, atol=1e-16)
    # qtgw = (qnew - q)/dt is a catastrophic cancellation: its noise
    # floor is ulp(q)/dt ~ 1.5e-25 per drag_prof call (~3e-25 for the
    # two-source sum); measured max abs diff 1.2e-25 at rel 1.1e-8.
    np.testing.assert_allclose(np.asarray(ptq), gold["chain_dqdt"],
                               rtol=RTOL, atol=1e-24)


# ---------------------------------------------------------------------------
# Tier-0 properties
# ---------------------------------------------------------------------------
def test_zero_heating_zero_tendency(gold, meta, params):
    """netdt = 0 everywhere => no convective waves, zero tendencies."""
    import jax.numpy as jnp
    bp = beres_params(gold, meta, params, True)
    src = gws.gw_beres_src(gold["lat"], gold["u"], gold["v"],
                           jnp.zeros_like(jnp.asarray(gold["netdt"])),
                           gold["zm"], params, bp)
    (src_l, tend_l, tau, ubm, ubi, xv, yv, c, hdepth, maxq0) = src
    assert float(jnp.abs(tau).max()) == 0.0
    assert float(jnp.abs(hdepth).max()) == 0.0
    out = run_drag(gold, meta, params,
                   (src_l, tend_l, tau, ubm, ubi, xv, yv, c),
                   meta["params"]["effgw_beres"], False)
    (_, utgw, vtgw, ttgw, qtgw, taucd, _, gwut, _, _) = out
    for a in (utgw, vtgw, ttgw, qtgw, taucd, gwut):
        assert float(jnp.abs(a).max()) == 0.0


def test_quiet_front_zero_source(gold, meta, params, fp):
    """frontgf below frontgfc => no frontal waves."""
    import jax.numpy as jnp
    quiet = jnp.zeros_like(jnp.asarray(gold["frontgf"]))
    _, _, tau, *_ = gws.gw_cm_src(gold["u"], gold["v"], quiet, params,
                                  fp)
    assert float(jnp.abs(tau).max()) == 0.0


def test_cm_source_spectral_symmetry(gold, meta, params):
    """The frontal source is symmetric in +/-l by construction (both
    Fortran golden and port), and wavenumber 0 is prohibited."""
    pgwv = meta["params"]["pgwv"]
    kbot = meta["params"]["kbotbg"]
    tau = gold["ctau"]
    np.testing.assert_array_equal(tau[:, :, kbot],
                                  tau[:, ::-1, kbot])
    assert np.all(tau[:, pgwv, :] == 0.0)


def test_zero_wind_spectral_symmetry(gold, meta, params):
    """Zero wind + symmetric source => tau stays symmetric in +/-l
    through the stress/tendency sweeps, gwut is antisymmetric, and the
    projected tendencies vanish (xv = yv = 0)."""
    import jax.numpy as jnp
    bp = beres_params(gold, meta, params, True)
    zero = jnp.zeros_like(jnp.asarray(gold["u"]))
    src = gws.gw_beres_src(gold["lat"], zero, zero, gold["netdt"],
                           gold["zm"], params, bp)
    (src_l, tend_l, tau, ubm, ubi, xv, yv, c, hdepth, _) = src
    assert float(jnp.abs(tau).max()) > 0.0     # something launched
    out = run_drag(gold, meta, params,
                   (src_l, tend_l, tau, ubm, ubi, xv, yv, c),
                   meta["params"]["effgw_beres"], False)
    (tau_o, utgw, vtgw, _, _, _, _, gwut, _, _) = out
    np.testing.assert_array_equal(np.asarray(tau_o),
                                  np.asarray(tau_o)[:, ::-1, :])
    np.testing.assert_array_equal(np.asarray(gwut),
                                  -np.asarray(gwut)[:, :, ::-1])
    assert float(jnp.abs(utgw).max()) == 0.0
    assert float(jnp.abs(vtgw).max()) == 0.0


def test_mec_momentum_balance(gold, meta, params):
    """momentum_energy_conservation adds exactly the stress remaining
    at tend_level to the column momentum budget."""
    g = meta["params"]["gravit"]
    dcol = (gold["mec_b_dudt"] - gold["bd_utgw"]) * gold["dpm"] / g
    tl = gold["btend_level_old"]
    taucd_tl = np.take_along_axis(
        gold["bd_taucd"], tl[:, None, None].repeat(4, axis=2),
        axis=1)[:, 0, :]
    expect = -(taucd_tl[:, 1] + taucd_tl[:, 0])   # east + west
    got = dcol.sum(axis=1)
    # measured residual 1.4e-17 vs stress scale 5e-2 (pure roundoff)
    assert (np.abs(got - expect)
            <= 1e-12 * np.abs(expect) + 1e-15).all()


def test_mec_energy_balance(gold, meta):
    """After the conservation fixer the column-integrated total energy
    tendency (s + KE terms, the quantity the routine zeroes) vanishes
    relative to its gross magnitude."""
    dt = meta["dt"]
    terms = gold["dpm"] * (
        gold["chain_dsdt"]
        + gold["chain_dudt"] * (gold["u"]
                                + gold["chain_dudt"] * 0.5 * dt)
        + gold["chain_dvdt"] * (gold["v"]
                                + gold["chain_dvdt"] * 0.5 * dt))
    resid = np.abs(terms.sum(axis=1))
    gross = np.abs(terms).sum(axis=1) + 1e-30
    assert (resid <= 1e-10 * gross + 1e-14).all()


def test_qtgw_column_mass_conservation(gold, meta):
    """The implicit GW diffusion has zero-flux boundaries (ksrf = 0,
    cc_top = 0), so the constituent tendency conserves column mass.
    Roundoff scale is the column tracer mass over dt (the solve
    perturbs q at the eps*q level); measured residual <= 2e-17 of it."""
    for tag in ("bd", "cd"):
        for m in range(gold["q"].shape[2]):
            col = (gold[f"{tag}_qtgw"][:, :, m] * gold["dpm"]).sum(axis=1)
            scale = ((np.abs(gold["q"][:, :, m]) * gold["dpm"])
                     .sum(axis=1) / meta["dt"])
            assert (np.abs(col) <= 1e-14 * scale).all()


def test_zero_source_tau_zero_tendency_prop(gold, meta, params):
    """tau = 0 through drag_prof => strictly zero wind tendency (the
    limiters never create stress)."""
    import jax.numpy as jnp
    inp = tuple(gold[f"c{k}"] for k in
                ("src_level", "tend_level", "tau", "ubm", "ubi", "xv",
                 "yv", "c"))
    src_l, tend_l, _, ubm, ubi, xv, yv, c = inp
    tau0 = jnp.zeros_like(jnp.asarray(gold["ctau"]))
    out = run_drag(gold, meta, params,
                   (src_l, tend_l, tau0, ubm, ubi, xv, yv, c),
                   meta["params"]["effgw_cm"], False)
    assert float(jnp.abs(out[1]).max()) == 0.0
    assert float(jnp.abs(out[4]).max()) == 0.0
