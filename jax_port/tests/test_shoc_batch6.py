"""Tier-0 tests for SHOC batch 6: third moments and the assumed PDF."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scream_jax.shoc import (  # noqa: E402
    clipping_diag_third_shoc_moments,
    constants as sc,
    diag_third_shoc_moments,
    shoc_assumed_pdf,
)
from scream_jax.shoc.assumed_pdf import (  # noqa: E402
    shoc_assumed_pdf_compute_s,
    shoc_assumed_pdf_vv_parameters,
)


def _grid(ncol=3, nlev=24, seed=0):
    rng = np.random.default_rng(seed)
    zi = np.linspace(12000.0, 0.0, nlev + 1)[None, :] * np.ones((ncol, 1))
    zt = 0.5 * (zi[:, :-1] + zi[:, 1:])
    dz_zt = zi[:, :-1] - zi[:, 1:]
    dz_zi = np.concatenate(
        [np.full((ncol, 1), 1.0), zt[:, :-1] - zt[:, 1:], zt[:, -1:]], axis=-1)
    return rng, zt, zi, dz_zt, dz_zi


def test_w3_clipping():
    w_sec_zi = np.full((2, 5), 0.5)
    bound = sc.w3clip * np.sqrt(2 * 0.5 ** 3)
    w3 = np.array([[0.0, bound * 0.99, bound * 1.01, -bound * 1.01, -bound * 0.99],
                   [0.1, -0.1, 3.0, -3.0, 0.0]])
    out = np.asarray(clipping_diag_third_shoc_moments(w_sec_zi, w3))
    # Values within the bound unchanged; outside replaced by 0.02
    np.testing.assert_allclose(out[0], [0.0, bound * 0.99, 0.02, 0.02, -bound * 0.99])
    np.testing.assert_allclose(out[1], [0.1, -0.1, 0.02, 0.02, 0.0])


def test_diag_third_moments_structure_and_symmetric_zero():
    rng, zt, zi, dz_zt, dz_zi = _grid()
    ncol, nlev = zt.shape
    dz_safe = np.where(dz_zi == 0, 1.0, dz_zi)

    # Uniform, symmetric turbulence with no gradients -> w3 = 0 everywhere
    w_sec = np.full((ncol, nlev), 0.3)
    tke = np.full((ncol, nlev), 0.45)
    thl_sec = np.full((ncol, nlev + 1), 0.1)
    wthl_sec = np.zeros((ncol, nlev + 1))
    iso = np.full((ncol, nlev), 100.0)
    brunt = np.zeros((ncol, nlev))
    thetal = np.full((ncol, nlev), 300.0)

    w3 = np.asarray(diag_third_shoc_moments(
        7.0, False, w_sec, thl_sec, wthl_sec, iso, brunt, thetal, tke,
        dz_zt, dz_safe, zt, zi))
    assert w3.shape == (ncol, nlev + 1)
    np.testing.assert_array_equal(w3[:, 0], 0.0)
    np.testing.assert_array_equal(w3[:, -1], 0.0)
    np.testing.assert_allclose(w3, 0.0, atol=1e-14)

    # 1.5-TKE closure: identically zero regardless of inputs
    w3b = np.asarray(diag_third_shoc_moments(
        7.0, True, w_sec, thl_sec, rng.uniform(-0.1, 0.1, (ncol, nlev + 1)),
        iso, brunt, thetal, tke, dz_zt, dz_safe, zt, zi))
    np.testing.assert_array_equal(w3b, 0.0)


def test_vv_parameters_moments_consistency():
    # The two-Gaussian mixture must reproduce zero mean and unit-ish
    # normalized variance: a*w1_1 + (1-a)*w1_2 = 0 (tilde), and
    # a*(w1_1^2+0.4) + (1-a)*(w1_2^2+0.4) = 1.
    rng = np.random.default_rng(1)
    w_sec = rng.uniform(0.05, 2.0, 50)
    w3 = rng.uniform(-1.0, 1.0, 50) * w_sec ** 1.5
    skew, w1_1, w1_2, w2_1, w2_2, a = (np.asarray(x) for x in
        shoc_assumed_pdf_vv_parameters(np.zeros(50), w_sec, w3))
    np.testing.assert_allclose(a * w1_1 + (1 - a) * w1_2, 0.0, atol=1e-12)
    np.testing.assert_allclose(a * (w1_1 ** 2) + (1 - a) * (w1_2 ** 2) + 0.4,
                               1.0, rtol=1e-12)
    assert np.all((a >= 0.01) & (a <= 0.99))


def test_compute_s_limits():
    # Saturated plume mean (qw1 >> qs): C -> 1; dry plume: C -> 0, qn = 0.
    p = np.full(4, 85000.0)
    qs = np.full(4, 8e-3)
    beta = np.full(4, 6e-6)
    qw1 = np.array([2e-2, 8e-3, 1e-3, 1e-4])   # very wet ... very dry
    thl2 = np.full(4, 0.04)
    qw2 = np.full(4, 1e-8)
    s, std_s, qn, C = (np.asarray(x) for x in shoc_assumed_pdf_compute_s(
        qw1, qs, beta, p, thl2, qw2, np.sqrt(thl2), np.sqrt(qw2), np.zeros(4)))
    assert C[0] > 0.99 and qn[0] > 0
    assert C[-1] < 0.01
    assert np.all((C >= 0) & (C <= 1)) and np.all(qn >= 0)
    # zero variance branch: s>0 -> C=1, qn=s; s<0 -> C=0, qn=0
    s2, _, qn2, C2 = (np.asarray(x) for x in shoc_assumed_pdf_compute_s(
        qw1, qs, beta, p, np.zeros(4), np.zeros(4), np.zeros(4), np.zeros(4),
        np.zeros(4)))
    pos = s2 > 0
    np.testing.assert_array_equal(C2[pos], 1.0)
    np.testing.assert_allclose(qn2[pos], s2[pos], rtol=1e-15)
    np.testing.assert_array_equal(C2[~pos], 0.0)
    np.testing.assert_array_equal(qn2[~pos], 0.0)


def test_assumed_pdf_driver_dry_and_moist():
    rng, zt, zi, dz_zt, dz_zi = _grid(seed=2)
    ncol, nlev = zt.shape
    nlevi = nlev + 1
    pres = np.linspace(300e2, 1000e2, nlev)[None] * np.ones((ncol, 1))

    common = dict(
        w_field=np.zeros((ncol, nlev)),
        thl_sec=np.full((ncol, nlevi), 0.09),
        qw_sec=np.full((ncol, nlevi), 1e-7),
        wthl_sec=np.full((ncol, nlevi), 0.01),
        w_sec=np.full((ncol, nlev), 0.2),
        wqw_sec=np.full((ncol, nlevi), 1e-5),
        qwthl_sec=np.full((ncol, nlevi), 1e-4),
        w3=np.full((ncol, nlevi), 0.05),
        pres=pres, zt_grid=zt, zi_grid=zi,
        shoc_ql_in=np.zeros((ncol, nlev)))

    # Dry mean state. NOTE: the prescribed flux moments imply plume
    # moisture excursions of ~1e-4 kg/kg, which DO saturate at the cold
    # upper levels (qs ~ 3e-5 at 213 K) — correct scheme behavior. The
    # no-cloud invariant only holds where it is warm (qs ~ 1e-2), i.e.
    # the lower half of the column.
    out_dry = shoc_assumed_pdf(
        np.full((ncol, nlev), 300.0), np.full((ncol, nlev), 1e-6),
        common["w_field"], common["thl_sec"], common["qw_sec"], 300.0, False,
        common["wthl_sec"], common["w_sec"], common["wqw_sec"],
        common["qwthl_sec"], common["w3"], pres, zt, zi, common["shoc_ql_in"])
    cldfrac, ql = np.asarray(out_dry[0]), np.asarray(out_dry[1])
    lower = slice(nlev // 2, None)
    np.testing.assert_allclose(cldfrac[:, lower], 0.0, atol=1e-10)
    np.testing.assert_allclose(ql[:, lower], 0.0, atol=1e-12)
    assert np.all((cldfrac >= 0) & (cldfrac <= 1))

    # Saturated lower levels: cloud must appear, with 0 <= cldfrac <= 1
    qw_moist = np.full((ncol, nlev), 1e-4)
    qw_moist[:, -6:] = 0.03
    thetal_cool = np.full((ncol, nlev), 285.0)
    out_moist = shoc_assumed_pdf(
        thetal_cool, qw_moist,
        common["w_field"], common["thl_sec"], common["qw_sec"], 300.0, True,
        common["wthl_sec"], common["w_sec"], common["wqw_sec"],
        common["qwthl_sec"], common["w3"], pres, zt, zi, common["shoc_ql_in"])
    (cldfrac2, ql2, wqls2, wthv2, sgs_ql2, cond2, evap2) = \
        (np.asarray(a) for a in out_moist)
    assert np.all((cldfrac2 >= 0) & (cldfrac2 <= 1))
    assert cldfrac2[:, -3:].min() > 0.5     # deep saturation -> cloudy
    assert ql2[:, -3:].min() > 1e-4
    assert np.all(ql2 >= 0) and np.all(sgs_ql2 >= 0)
    # cond/evap diagnostics: with ql_in = 0, evap = 0 and cond = ql/dt
    np.testing.assert_allclose(cond2, ql2 / 300.0, rtol=1e-12)
    np.testing.assert_array_equal(evap2, 0.0)
    assert np.all(np.isfinite(wqls2)) and np.all(np.isfinite(wthv2))
