"""Tier-0 tests for SHOC batch 3: length-scale and TKE chains.

Invariants mirror the C++ run_property tests (shoc_length_tests.cpp,
shoc_tke_tests.cpp, shoc_eddy_diffusivities_tests.cpp, ...).
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scream_jax.foundation import constants as c  # noqa: E402
from scream_jax.shoc import (  # noqa: E402
    adv_sgs_tke,
    compute_brunt_shoc_length,
    compute_l_inf_shoc_length,
    compute_shr_prod,
    constants as sc,
    eddy_diffusivities,
    integ_column_stability,
    isotropic_ts,
    shoc_length,
    shoc_tke,
)


def _grid(ncol=3, nlev=32, seed=0):
    rng = np.random.default_rng(seed)
    zi = np.linspace(20000.0, 0.0, nlev + 1)[None, :] * np.ones((ncol, 1))
    zt = 0.5 * (zi[:, :-1] + zi[:, 1:])
    dz_zt = zi[:, :-1] - zi[:, 1:]
    dz_zi = np.concatenate(
        [np.zeros((ncol, 1)), zt[:, :-1] - zt[:, 1:], zt[:, -1:]], axis=-1)
    return rng, zt, zi, dz_zt, dz_zi


def test_brunt_sign_from_stratification():
    _, zt, zi, dz_zt, _ = _grid()
    # thv increasing with height (stable) -> brunt > 0
    thv = 300.0 + 0.005 * zt
    thv_zi = 300.0 + 0.005 * zi
    brunt = np.asarray(compute_brunt_shoc_length(dz_zt, thv, thv_zi))
    assert np.all(brunt > 0)
    # unstable
    brunt2 = np.asarray(compute_brunt_shoc_length(dz_zt, 300.0 - 0.005 * zt,
                                                  300.0 - 0.005 * zi))
    assert np.all(brunt2 < 0)


def test_l_inf_is_weighted_height():
    _, zt, _, dz_zt, _ = _grid(seed=1)
    tke = np.full(zt.shape, 0.25)
    # Uniform tke -> l_inf = 0.1 * dz-weighted mean height
    l_inf = np.asarray(compute_l_inf_shoc_length(zt, dz_zt, tke))
    expected = 0.1 * (zt * dz_zt).sum(-1) / dz_zt.sum(-1)
    np.testing.assert_allclose(l_inf, expected, rtol=1e-13)


def test_shoc_length_bounds_and_stability_shrink():
    rng, zt, zi, dz_zt, _ = _grid(seed=2)
    tke = rng.uniform(sc.mintke, 1.0, zt.shape)
    dx = dy = np.full(zt.shape[0], 5000.0)

    thv_stable = 300.0 + 0.01 * zt
    _, mix_stable = (np.asarray(a) for a in shoc_length(
        0.5, dx, dy, zt, zi, dz_zt, tke, thv_stable))
    thv_neutral = np.full(zt.shape, 300.0)
    _, mix_neutral = (np.asarray(a) for a in shoc_length(
        0.5, dx, dy, zt, zi, dz_zt, tke, thv_neutral))

    for mix in (mix_stable, mix_neutral):
        assert np.all(mix >= sc.minlen) and np.all(mix <= 5000.0)
    # Stable stratification must not increase the mixing length anywhere
    assert np.all(mix_stable <= mix_neutral + 1e-12)


def test_integ_column_stability_pressure_gate():
    _, zt, _, dz_zt, _ = _grid(seed=3)
    brunt = np.full(zt.shape, 1e-4)
    # All pressure above 800 hPa threshold -> full integral; below -> zero
    p_hi = np.full(zt.shape, 90000.0)
    p_lo = np.full(zt.shape, 50000.0)
    full = np.asarray(integ_column_stability(dz_zt, p_hi, brunt))
    np.testing.assert_allclose(full, (dz_zt * 1e-4).sum(-1), rtol=1e-13)
    np.testing.assert_array_equal(
        np.asarray(integ_column_stability(dz_zt, p_lo, brunt)), 0.0)


def test_shr_prod_boundaries_and_positivity():
    rng, zt, _, _, dz_zi = _grid(seed=4)
    u = rng.uniform(-30, 30, zt.shape)
    v = rng.uniform(-30, 30, zt.shape)
    dz_safe = np.where(dz_zi == 0, 1.0, dz_zi)
    sterm = np.asarray(compute_shr_prod(dz_safe, u, v))
    assert sterm.shape[-1] == zt.shape[-1] + 1
    np.testing.assert_array_equal(sterm[:, 0], 0.0)
    np.testing.assert_array_equal(sterm[:, -1], 0.0)
    assert np.all(sterm >= 0)
    # No shear -> no production
    sterm0 = np.asarray(compute_shr_prod(dz_safe, np.full(zt.shape, 7.0),
                                         np.full(zt.shape, -3.0)))
    np.testing.assert_array_equal(sterm0, 0.0)


def test_adv_sgs_tke_bounds_and_dissipation():
    rng, zt, _, _, _ = _grid(seed=5)
    shape = zt.shape
    mix = rng.uniform(sc.minlen, 1000.0, shape)
    tke0 = rng.uniform(sc.mintke, 2.0, shape)
    tk = rng.uniform(0.0, 10.0, shape)
    sterm = rng.uniform(0.0, 1e-3, shape)
    wthv = rng.uniform(-0.05, 0.05, shape)
    brunt = rng.uniform(-1e-4, 1e-4, shape)

    tke1, a_diss = (np.asarray(a) for a in adv_sgs_tke(
        300.0, False, mix, wthv, sterm, tk, brunt, tke0))
    assert np.all(tke1 >= sc.mintke) and np.all(tke1 <= sc.maxtke)
    assert np.all(a_diss >= 0)
    # Pure dissipation (no production): tke must not grow
    tke2, _ = (np.asarray(a) for a in adv_sgs_tke(
        300.0, False, mix, np.zeros(shape), np.zeros(shape),
        np.zeros(shape), brunt, tke0))
    assert np.all(tke2 <= tke0 + 1e-15)
    # 1.5-TKE closure: unstable brunt (<0) acts as production
    tke3, _ = (np.asarray(a) for a in adv_sgs_tke(
        300.0, True, mix, np.zeros(shape), np.zeros(shape),
        np.zeros(shape), np.full(shape, -1e-3), tke0))
    assert np.all(tke3 >= tke2 - 1e-15)


def test_isotropic_ts_caps_and_stability_damping():
    shape = (2, 8)
    tke = np.full(shape, 0.4)
    a_diss = np.full(shape, 1e-6)      # huge tscale -> hits maxiso when lam=0
    brunt_neg = np.full(shape, -1e-4)  # lambda zeroed
    iso = np.asarray(isotropic_ts(0.001, 0.08, 2.65, 0.02,
                                  np.zeros(2), tke, a_diss, brunt_neg))
    np.testing.assert_array_equal(iso, sc.maxiso)
    # Stable brunt damps the timescale below the undamped value
    brunt_pos = np.full(shape, 1e-4)
    iso2 = np.asarray(isotropic_ts(0.001, 0.08, 2.65, 0.02,
                                   np.full(2, 10.0), tke, np.full(shape, 1e-3),
                                   brunt_pos))
    undamped = 2 * 0.4 / 1e-3
    assert np.all(iso2 < undamped)


def test_eddy_diffusivities_regimes():
    rng, zt, _, _, _ = _grid(seed=6)
    shape = zt.shape
    mix = rng.uniform(sc.minlen, 1000.0, shape)
    sterm = rng.uniform(1e-6, 1e-3, shape)
    iso = rng.uniform(1.0, 1000.0, shape)
    tke = rng.uniform(sc.mintke, 2.0, shape)
    pblh = np.full(shape[0], 1000.0)

    # Normal regime: tkh/tk = Ckh/Ckm * isotropy * tke everywhere
    tabs_warm = np.full(shape, 280.0)
    tkh, tk = (np.asarray(a) for a in eddy_diffusivities(
        0.1, 0.1, pblh, zt, tabs_warm, mix, sterm, iso, tke))
    np.testing.assert_allclose(tkh, 0.1 * iso * tke, rtol=1e-14)
    np.testing.assert_allclose(tk, 0.1 * iso * tke, rtol=1e-14)

    # Runaway-cooling regime: below pblh+200 use stable-PBL form
    tabs_cold = np.full(shape, 280.0)
    tabs_cold[:, -1] = 150.0
    tkh2, _ = (np.asarray(a) for a in eddy_diffusivities(
        0.1, 0.1, pblh, zt, tabs_cold, mix, sterm, iso, tke))
    low = zt < (pblh[:, None] + 200.0)
    np.testing.assert_allclose(tkh2[low], (0.1 * mix ** 2 * np.sqrt(sterm))[low],
                               rtol=1e-14)
    np.testing.assert_allclose(tkh2[~low], (0.1 * iso * tke)[~low], rtol=1e-14)


def test_shoc_tke_driver_runs_and_bounds():
    rng, zt, zi, dz_zt, dz_zi = _grid(seed=7)
    shape = zt.shape
    dz_safe = np.where(dz_zi == 0, 1.0, dz_zi)
    tke, tk, tkh, iso = (np.asarray(a) for a in shoc_tke(
        300.0, 0.001, 0.08, 2.65, 0.02, 0.1, 0.1, False,
        rng.uniform(-0.05, 0.05, shape),          # wthv_sec
        rng.uniform(sc.minlen, 1000.0, shape),    # shoc_mix
        dz_safe, dz_zt,
        np.linspace(300e2, 1000e2, shape[1])[None] * np.ones((shape[0], 1)),
        np.full(shape, 280.0),                    # tabs
        rng.uniform(-20, 20, shape), rng.uniform(-20, 20, shape),
        rng.uniform(-1e-4, 1e-4, shape),          # brunt
        zt, zi, np.full(shape[0], 1000.0),
        rng.uniform(sc.mintke, 1.0, shape),       # tke
        rng.uniform(0.0, 10.0, shape)))           # tk (previous)
    assert np.all((tke >= sc.mintke) & (tke <= sc.maxtke))
    assert np.all(iso <= sc.maxiso) and np.all(iso >= 0)
    assert np.all(np.isfinite(tkh)) and np.all(np.isfinite(tk))
    assert np.all(tkh >= 0) and np.all(tk >= 0)
