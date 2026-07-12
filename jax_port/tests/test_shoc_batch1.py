"""Tier-0 tests for SHOC batch 1: grid, thermo, interp, check_tke.

Property tests mirror the C++ run_property tests in
components/eamxx/src/physics/shoc/tests/ (shoc_grid_tests.cpp,
shoc_linear_interp_tests.cpp, shoc_compute_shoc_temperature_tests.cpp, ...);
golden whole-scheme validation happens at the shoc_main level.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scream_jax.foundation import constants as c  # noqa: E402
from scream_jax.shoc import (  # noqa: E402
    check_tke,
    compute_shoc_temperature,
    compute_shoc_vapor,
    compute_tmpi,
    constants as sc,
    dp_inverse,
    linear_interp,
    shoc_grid,
)


def _column_grid(ncol=4, nlev=72, seed=0):
    rng = np.random.default_rng(seed)
    # Monotonically decreasing interface heights, surface at z=0.
    zi = np.sort(rng.uniform(30.0, 40000.0, (ncol, nlev - 1)), axis=-1)[:, ::-1]
    zi = np.concatenate([np.full((ncol, 1), 42000.0), zi, np.zeros((ncol, 1))], axis=-1)
    zt = 0.5 * (zi[:, :-1] + zi[:, 1:])
    pdel = rng.uniform(100.0, 2500.0, (ncol, nlev))
    return zt, zi, pdel


def test_shoc_grid_properties():
    # Mirrors shoc_grid_tests run_property: dz > 0, rho physical, and
    # column sums telescope back to the column height.
    zt, zi, pdel = _column_grid()
    dz_zt, dz_zi, rho_zt = (np.asarray(a) for a in shoc_grid(zt, zi, pdel))

    assert np.all(dz_zt > 0)
    np.testing.assert_allclose(dz_zt.sum(axis=-1), zi[:, 0], rtol=1e-12)

    # dz_zi boundaries: 0 at top, zt(nlev-1) at surface; interior positive.
    assert np.all(dz_zi[:, 0] == 0.0)
    np.testing.assert_array_equal(dz_zi[:, -1], zt[:, -1])
    assert np.all(dz_zi[:, 1:-1] > 0)
    np.testing.assert_allclose(dz_zi[:, 1:-1], zt[:, :-1] - zt[:, 1:], rtol=0)

    # Density: rho = pdel/(g dz), all positive and < 2 kg/m3 for sane inputs.
    np.testing.assert_allclose(rho_zt, pdel / (c.gravit * dz_zt), rtol=1e-15)
    assert np.all(rho_zt > 0)


def test_dp_inverse_and_tmpi():
    zt, zi, pdel = _column_grid(seed=1)
    dz_zt, dz_zi, rho_zt = shoc_grid(zt, zi, pdel)
    rdp = np.asarray(dp_inverse(rho_zt, dz_zt))
    # 1/(g rho dz) = 1/pdel by construction
    np.testing.assert_allclose(rdp, 1.0 / pdel, rtol=1e-13)

    rho_zi = np.random.default_rng(2).uniform(0.1, 1.4, dz_zi.shape)
    dz_zi_safe = np.where(np.asarray(dz_zi) == 0, 1.0, np.asarray(dz_zi))
    tmpi = np.asarray(compute_tmpi(300.0, rho_zi, dz_zi_safe))
    assert np.all(tmpi[:, 0] == 0.0)
    np.testing.assert_allclose(
        tmpi[:, 1:], 300.0 * c.gravit * rho_zi[:, 1:] / dz_zi_safe[:, 1:], rtol=1e-15)


def test_vapor_and_temperature_roundtrip():
    rng = np.random.default_rng(3)
    shape = (5, 72)
    ql = rng.uniform(0.0, 2e-3, shape)
    qw = ql + rng.uniform(1e-4, 2e-2, shape)
    qv = np.asarray(compute_shoc_vapor(qw, ql))
    np.testing.assert_array_equal(qv, qw - ql)
    assert np.all(qv > 0)

    # compute_shoc_temperature: with ql=0, tabs = thetal/inv_exner; adding
    # ql adds exactly (Lv/cp) ql (mirrors run_property points in
    # shoc_compute_shoc_temperature_tests.cpp).
    thetal = rng.uniform(250.0, 320.0, shape)
    inv_exner = rng.uniform(1.0, 2.5, shape)
    t0 = np.asarray(compute_shoc_temperature(thetal, np.zeros(shape), inv_exner))
    np.testing.assert_allclose(t0, thetal / inv_exner, rtol=1e-15)
    t1 = np.asarray(compute_shoc_temperature(thetal, ql, inv_exner))
    # atol accounts for cancellation: t1-t0 isolates a ~mK term from ~300 K.
    np.testing.assert_allclose(t1 - t0, (c.LatVap / c.CP) * ql,
                               rtol=1e-9, atol=1e-10)


def test_linear_interp_exact_on_linear_functions():
    # Both directions must reproduce a linear function exactly (up to fp),
    # including the extrapolating end points of the zt->zi direction.
    nlev = 12
    rng = np.random.default_rng(4)
    zi = np.sort(rng.uniform(0.0, 1000.0, nlev + 1))[::-1].copy()
    zt = 0.5 * (zi[:-1] + zi[1:])
    a, b = 3.7, -42.0

    # zi -> zt (km1 = nlev+1 = km2+1)
    y_zt = np.asarray(linear_interp(zi, zt, a * zi + b, nlev + 1, nlev, sc.largeneg))
    np.testing.assert_allclose(y_zt, a * zt + b, rtol=1e-12)

    # zt -> zi (km2 = nlev+1 = km1+1), extrapolation at both ends
    y_zi = np.asarray(linear_interp(zt, zi, a * zt + b, nlev, nlev + 1, sc.largeneg))
    np.testing.assert_allclose(y_zi, a * zi + b, rtol=1e-12)


def test_linear_interp_minthresh_floor():
    x1 = np.array([3.0, 2.0, 1.0, 0.0])
    y1 = np.array([-5.0, -4.0, -3.0, -2.0])
    x2 = 0.5 * (x1[:-1] + x1[1:])
    y2 = np.asarray(linear_interp(x1, x2, y1, 4, 3, -3.5))
    assert np.all(y2 >= -3.5)
    # unclipped values would be [-4.5, -3.5, -2.5]
    np.testing.assert_allclose(y2, [-3.5, -3.5, -2.5], rtol=1e-14)


def test_check_tke_clips_only_below_min():
    tke = np.array([[0.0, -1.0, sc.mintke, 2 * sc.mintke, 5.0]])
    out = np.asarray(check_tke(tke))
    np.testing.assert_array_equal(out, [[sc.mintke, sc.mintke, sc.mintke,
                                         2 * sc.mintke, 5.0]])
