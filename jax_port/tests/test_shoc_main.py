"""Tier-0 smoke/invariant tests for the assembled shoc_main.

The authoritative validation is the Tier-1 comparison against
jax_port/golden/shoc_218x72_dt1800_5steps.npz, which additionally requires
the process-interface pre/post conversions (see STATUS.md).
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scream_jax.foundation import constants as c  # noqa: E402
from scream_jax.shoc import constants as sc, shoc_init, shoc_main  # noqa: E402


def test_shoc_init_npbl():
    pref = np.linspace(2000.0, 101000.0, 72)
    npbl = shoc_init(72, 0, pref)
    expected = 72 - np.argmax(pref >= sc.pblmaxp)
    assert npbl == expected
    # All levels above pblmaxp -> npbl = 1
    assert shoc_init(72, 0, np.full(72, 1e3)) == 1


def _state(ncol=4, nlev=48, nq=3, seed=0):
    rng = np.random.default_rng(seed)
    zi = np.linspace(28000.0, 0.0, nlev + 1)[None, :] * np.ones((ncol, 1))
    zt = 0.5 * (zi[:, :-1] + zi[:, 1:])
    # Hydrostatic-ish pressure profile
    p0 = 101325.0
    presi = p0 * np.exp(-zi / 8000.0)
    pres = p0 * np.exp(-zt / 8000.0)
    pdel = presi[:, 1:] - presi[:, :-1]
    inv_exner = (p0 / pres) ** (c.Rair / c.CP)

    thetal = 290.0 + 30.0 * (zt / zt[:, :1]) + rng.uniform(-0.5, 0.5, zt.shape)
    qw = np.clip(0.016 * np.exp(-zt / 3000.0), 1e-7, None)
    ql = np.zeros_like(qw)
    tke = np.full(zt.shape, 0.01)
    u = 5.0 + rng.uniform(-1, 1, zt.shape)
    v = -3.0 + rng.uniform(-1, 1, zt.shape)
    host_dse = c.CP * (thetal / inv_exner) + c.gravit * zt
    return dict(
        zt=zt, zi=zi, pres=pres, presi=presi, pdel=pdel, inv_exner=inv_exner,
        thetal=thetal, qw=qw, ql=ql, tke=tke, u=u, v=v, host_dse=host_dse,
        thv=thetal * (1 + 0.61 * qw),
        w_field=np.zeros_like(zt),
        qtr=rng.uniform(0.0, 1e-3, (ncol, nq, nlev)),
        wthv_sec=np.full(zt.shape, 1e-3),
        tk=np.full(zt.shape, 1.0),
        cldfrac=np.zeros_like(zt),
        ncol=ncol, nlev=nlev, nq=nq)


def test_shoc_main_runs_and_respects_invariants():
    s = _state()
    ncol, nlev, nq = s["ncol"], s["nlev"], s["nq"]
    zeros = np.zeros(ncol)
    npbl = shoc_init(nlev, 0, np.asarray(s["pres"]).mean(axis=0))

    out = shoc_main(
        300.0, 2, npbl,
        0.001, 0.08, 2.65, 0.02,      # lambda_*
        1.0, 1.0, 1.0, 1.0,           # tune factors
        0.5, 7.0, 0.1, 0.1,           # length_fac, c_diag_3rd_mom, Ckh, Ckm
        False, True,                  # shoc_1p5tke, extra_diags
        np.full(ncol, 10000.0), np.full(ncol, 10000.0),
        s["zt"], s["zi"], s["pres"], s["presi"], s["pdel"], s["thv"],
        s["w_field"],
        np.full(ncol, 0.03), np.full(ncol, 1e-5),     # wthl, wqw
        np.full(ncol, 0.05), np.full(ncol, -0.02),    # uw, vw
        np.zeros((ncol, nq)), s["inv_exner"], zeros,
        s["host_dse"], s["tke"], s["thetal"], s["qw"], s["u"], s["v"],
        s["wthv_sec"], s["qtr"], s["tk"], s["cldfrac"], s["ql"])

    for k, v in out.items():
        assert np.all(np.isfinite(np.asarray(v))), f"non-finite in {k}"

    tke = np.asarray(out["tke"])
    assert np.all((tke >= sc.mintke) & (tke <= sc.maxtke))
    assert np.all(np.asarray(out["shoc_ql"]) >= 0)
    cf = np.asarray(out["shoc_cldfrac"])
    assert np.all((cf >= 0) & (cf <= 1))
    mix = np.asarray(out["shoc_mix"])
    assert np.all((mix >= sc.minlen) & (mix <= 10000.0))
    pblh = np.asarray(out["pblh"])
    assert np.all((pblh > 0) & (pblh < 28000.0))
    assert np.all(np.asarray(out["ustar"]) >= sc.ustar_min)
    # Boundary structure of second moments survived the loop
    np.testing.assert_array_equal(np.asarray(out["wthl_sec"])[:, 0], 0.0)
    np.testing.assert_allclose(np.asarray(out["wthl_sec"])[:, -1], 0.03,
                               rtol=1e-13)

    # Column water (qw) is conserved by the internal loop except for the
    # surface flux source: gain = nadv * dt * wqw * g * rho_zi_sfc
    from scream_jax.shoc import linear_interp, shoc_grid
    dz_zt, _, rho_zt = (np.asarray(a) for a in
                        shoc_grid(s["zt"], s["zi"], s["pdel"]))
    rho_zi = np.asarray(linear_interp(s["zt"], s["zi"], rho_zt,
                                      nlev, nlev + 1, 0.0))
    w = s["pdel"] / c.gravit
    gain = (np.asarray(out["qw"]) * w).sum(-1) - (s["qw"] * w).sum(-1)
    expected = 2 * 300.0 * 1e-5 * rho_zi[:, -1]
    np.testing.assert_allclose(gain, expected, rtol=1e-9)


def test_shoc_main_zero_flux_conserves_energy_integral():
    s = _state(seed=1)
    ncol, nlev, nq = s["ncol"], s["nlev"], s["nq"]
    zeros = np.zeros(ncol)
    npbl = shoc_init(nlev, 0, np.asarray(s["pres"]).mean(axis=0))

    out = shoc_main(
        300.0, 1, npbl,
        0.001, 0.08, 2.65, 0.02, 1.0, 1.0, 1.0, 1.0, 0.5, 7.0, 0.1, 0.1,
        False, False,
        np.full(ncol, 10000.0), np.full(ncol, 10000.0),
        s["zt"], s["zi"], s["pres"], s["presi"], s["pdel"], s["thv"],
        s["w_field"], zeros, zeros, zeros, zeros,
        np.zeros((ncol, nq)), s["inv_exner"], zeros,
        s["host_dse"], s["tke"], s["thetal"], s["qw"], s["u"], s["v"],
        s["wthv_sec"], s["qtr"], s["tk"], s["cldfrac"], s["ql"])

    # With zero surface fluxes the energy fixer restores the total
    # (dse + ke + L-weighted water) integral to its pre-SHOC value.
    w = s["pdel"] / c.gravit
    def total(dse, qw, ql, u, v):
        se = (dse * w).sum(-1)
        ke = (0.5 * (u ** 2 + v ** 2) * w).sum(-1)
        wv = ((qw - ql) * w).sum(-1)
        wl = (ql * w).sum(-1)
        return se + ke + (c.LatVap + c.LatIce) * wv + c.LatIce * wl

    before = total(s["host_dse"], s["qw"], s["ql"], s["u"], s["v"])
    after = total(np.asarray(out["host_dse"]), np.asarray(out["qw"]),
                  np.asarray(out["shoc_ql"]), np.asarray(out["u_wind"]),
                  np.asarray(out["v_wind"]))
    np.testing.assert_allclose(after, before, rtol=1e-10)
