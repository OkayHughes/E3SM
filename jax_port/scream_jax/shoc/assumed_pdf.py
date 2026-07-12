"""Assumed double-Gaussian PDF — the centerpiece of SHOC.

Computes SGS cloud fraction, liquid water, liquid-water flux, and buoyancy
flux following the appendix of Larson et al. (2002), Analytic Double
Gaussian 1.

Sources (components/eamxx/src/physics/shoc/impl/):
  shoc_assumed_pdf_impl.hpp (driver) and the 12 shoc_assumed_pdf_*_impl.hpp
  helper kernels.

Deviations from the C++ (all value-preserving for defined inputs):
- Pack masked-set becomes jnp.where; "compute the second plume only if it
  differs" shortcuts are evaluated unconditionally and selected (equal
  inputs give equal outputs).
- In compute_s the C++ leaves C/qn uninitialized when std_s is tiny and
  s <= 0; this port zero-initializes them (the value any non-positive
  garbage would be clamped to by the final qn<=0 clamp).
- Divisions the C++ performs only under a mask use a safe denominator on
  masked-out lanes (selected away); masked-in arithmetic is unchanged.
- NaN/min-temperature guards keep the clamps but not the printfs.
"""

import functools

import jax
import jax.numpy as jnp

from ..foundation import constants as c
from ..foundation.saturation import murphy_koop_svp
from . import constants as sc
from .interp import linear_interp

# Tolerances/thresholds from the driver (shoc_assumed_pdf_impl.hpp)
_THL_TOL = 1e-2
_RT_TOL = 1e-4
_W_TOL_SQD = 4e-4
_W_THRESH = 0.0
_TL_MIN = 100.0


def shoc_assumed_pdf_vv_parameters(w_first, w_sec, w3var):
    """Vertical-velocity PDF parameters
    (Functions::shoc_assumed_pdf_vv_parameters). Returns
    (Skew_w, w1_1, w1_2, w2_1, w2_2, a); w1_* still in tilde form."""
    cond = w_sec > _W_TOL_SQD
    tmp_val = 0.4
    one_m = 1.0 - tmp_val
    sqrtw2t = jnp.sqrt(1.0 - tmp_val)

    w_sec_safe = jnp.where(cond, w_sec, 1.0)
    skew_w = jnp.where(cond, w3var / jnp.sqrt(w_sec_safe ** 3), 0.0)
    a = jnp.where(
        cond,
        jnp.maximum(0.01, jnp.minimum(
            0.99,
            0.5 * (1.0 - skew_w * jnp.sqrt(
                1.0 / (4.0 * one_m ** 3 + skew_w ** 2))))),
        0.5)
    w1_1 = jnp.where(cond, jnp.sqrt((1.0 - a) / a) * sqrtw2t, w_first)
    w1_2 = jnp.where(cond, -jnp.sqrt(a / (1.0 - a)) * sqrtw2t, w_first)
    w2_1 = jnp.where(cond, tmp_val * w_sec, 0.0)
    w2_2 = jnp.where(cond, tmp_val * w_sec, 0.0)
    return skew_w, w1_1, w1_2, w2_1, w2_2, a


def _pdf_scalar_parameters(wxsec, sqrtw2, sqrtx, xsec, x_first,
                           w1_1, w1_2, skew_x, a, cond, default_first):
    """Shared algebra of the thl/qw parameter kernels (their masked-set
    bodies are identical; only the skewness input and tolerances differ).

    Returns (x1_1, x1_2, x2_1, x2_2, sqrtx2_1, sqrtx2_2).
    """
    corrtest = jnp.maximum(-1.0, jnp.minimum(1.0, wxsec / (sqrtw2 * sqrtx)))
    w1_1_safe = jnp.where(cond, w1_1, 1.0)
    w1_2_safe = jnp.where(cond, w1_2, 1.0)
    tmp1 = -corrtest / w1_1_safe   # tmp_val_1
    tmp2 = -corrtest / w1_2_safe   # tmp_val_2

    common = 1.0 - a * tmp2 ** 2 - (1.0 - a) * tmp1 ** 2
    skew_term = skew_x - a * tmp2 ** 3 - (1.0 - a) * tmp1 ** 3
    denom = jnp.where(cond, tmp1 - tmp2, 1.0)

    x2_1 = jnp.where(
        cond,
        jnp.minimum(100.0, jnp.maximum(
            0.0, (3.0 * tmp1 * common - skew_term) / (3.0 * a * denom))) * xsec,
        0.0)
    x2_2 = jnp.where(
        cond,
        jnp.minimum(100.0, jnp.maximum(
            0.0, (-3.0 * tmp2 * common + skew_term)
            / (3.0 * (1.0 - a) * denom))) * xsec,
        0.0)

    x1_1 = jnp.where(cond, tmp2 * sqrtx + x_first, default_first)
    x1_2 = jnp.where(cond, tmp1 * sqrtx + x_first, default_first)

    sqrtx2_1 = jnp.where(cond, jnp.sqrt(x2_1), 0.0)
    sqrtx2_2 = jnp.where(cond, jnp.sqrt(x2_2), 0.0)
    return x1_1, x1_2, x2_1, x2_2, sqrtx2_1, sqrtx2_2


def shoc_assumed_pdf_thl_parameters(wthlsec, sqrtw2, sqrtthl, thlsec,
                                    thl_first, w1_1, w1_2, skew_w, a):
    """thetal PDF parameters (Functions::shoc_assumed_pdf_thl_parameters).

    Skew_thl is 0 unless the compile-time dothetal_skew switch is on
    (shoc_constants.hpp has it False).
    """
    cond = (thlsec > _THL_TOL ** 2) & (jnp.abs(w1_2 - w1_1) > _W_THRESH)

    if sc.dothetal_skew:
        corrtest = jnp.maximum(-1.0, jnp.minimum(1.0, wthlsec / (sqrtw2 * sqrtthl)))
        w1_1_safe = jnp.where(cond, w1_1, 1.0)
        w1_2_safe = jnp.where(cond, w1_2, 1.0)
        tsign = jnp.abs(-corrtest / w1_1_safe - (-corrtest / w1_2_safe))
        skew_thl = jnp.where(
            tsign > 0.4, 1.2 * skew_w,
            jnp.where(tsign > 0.2, ((1.2 * skew_w) / 0.2) * (tsign - 0.2), 0.0))
    else:
        skew_thl = jnp.zeros_like(skew_w)

    return _pdf_scalar_parameters(wthlsec, sqrtw2, sqrtthl, thlsec, thl_first,
                                  w1_1, w1_2, skew_thl, a, cond, thl_first)


def shoc_assumed_pdf_qw_parameters(wqwsec, sqrtw2, skew_w, sqrtqt, qwsec,
                                   w1_2, w1_1, qw_first, a):
    """Total-water PDF parameters (Functions::shoc_assumed_pdf_qw_parameters).

    Unlike thl, the qw skewness ramp is always active.
    """
    cond = (qwsec > _RT_TOL ** 2) & (jnp.abs(w1_2 - w1_1) > _W_THRESH)

    corrtest = jnp.maximum(-1.0, jnp.minimum(1.0, wqwsec / (sqrtw2 * sqrtqt)))
    w1_1_safe = jnp.where(cond, w1_1, 1.0)
    w1_2_safe = jnp.where(cond, w1_2, 1.0)
    tsign = jnp.abs(-corrtest / w1_1_safe - (-corrtest / w1_2_safe))
    skew_qw = jnp.where(
        tsign > 0.4, 1.2 * skew_w,
        jnp.where(tsign > 0.2, ((1.2 * skew_w) / 0.2) * (tsign - 0.2), 0.0))

    return _pdf_scalar_parameters(wqwsec, sqrtw2, sqrtqt, qwsec, qw_first,
                                  w1_1, w1_2, skew_qw, a, cond, qw_first)


def shoc_assumed_pdf_tilde_to_real(w_first, sqrtw2, w1):
    """Normalized (tilde) plume velocity -> physical: w1*sqrt(w2) + w_first."""
    return w1 * sqrtw2 + w_first


def shoc_assumed_pdf_inplume_correlations(sqrtqw2_1, sqrtthl2_1, a, sqrtqw2_2,
                                          sqrtthl2_2, qwthlsec, qw1_1, qw_first,
                                          thl1_1, thl_first, qw1_2, thl1_2):
    """Within-plume qw-thl correlation
    (Functions::shoc_assumed_pdf_inplume_correlations)."""
    testvar = a * sqrtqw2_1 * sqrtthl2_1 + (1.0 - a) * sqrtqw2_2 * sqrtthl2_2
    nonzero = testvar != 0.0
    testvar_safe = jnp.where(nonzero, testvar, 1.0)
    r = (qwthlsec - a * (qw1_1 - qw_first) * (thl1_1 - thl_first)
         - (1.0 - a) * (qw1_2 - qw_first) * (thl1_2 - thl_first)) / testvar_safe
    return jnp.where(nonzero, jnp.maximum(-1.0, jnp.minimum(1.0, r)), 0.0)


def shoc_assumed_pdf_compute_temperature(thl1, pval):
    """Plume liquid-water temperature: Tl1 = thl1 / (P0/p)^(Rd/cp)."""
    return thl1 / ((c.P0 / pval) ** (c.Rair / c.CP))


def shoc_assumed_pdf_compute_qs(Tl1_1, Tl1_2, pval):
    """Saturation mixing ratio and beta for both plumes
    (Functions::shoc_assumed_pdf_compute_qs; Murphy-Koop over liquid).

    Returns (qs1, beta1, qs2, beta2). The second plume is recomputed only
    where Tl1_1 != Tl1_2 in the C++; unconditional evaluation is identical.
    """
    def one(tl):
        es = murphy_koop_svp(tl, ice=False)
        qs = 0.622 * es / jnp.maximum(es, pval - es)
        beta = (c.Rair / c.RV) * (c.LatVap / (c.Rair * tl)) * (c.LatVap / (c.CP * tl))
        return qs, beta

    qs1, beta1 = one(Tl1_1)
    qs2, beta2 = one(Tl1_2)
    return qs1, beta1, qs2, beta2


def shoc_assumed_pdf_compute_s(qw1, qs, beta, pval, thl2, qw2,
                               sqrtthl2, sqrtqw2, r_qwthl):
    """Saturation-excess statistics for one plume
    (Functions::shoc_assumed_pdf_compute_s).

    Returns (s, std_s, qn, C): mean excess, its std dev, in-plume
    condensate, and cloud fraction.
    """
    sqrt2 = jnp.sqrt(2.0)
    sqrt2pi = jnp.sqrt(2.0 * jnp.pi)

    cthl = ((1.0 + beta * qw1) / (1.0 + beta * qs) ** 2) * (c.CP / c.LatVap) \
        * beta * qs * (pval / c.P0) ** (c.Rair / c.CP)
    cqt = 1.0 / (1.0 + beta * qs)

    std_s = jnp.sqrt(jnp.maximum(
        0.0,
        cthl ** 2 * thl2 + cqt ** 2 * qw2
        - 2.0 * cthl * sqrtthl2 * cqt * sqrtqw2 * r_qwthl))
    s = qw1 - qs * ((1.0 + beta * qw1) / (1.0 + beta * qs))

    tiny_floor = jnp.sqrt(jnp.finfo(jnp.float64).tiny) * 100.0
    not_small = std_s > tiny_floor
    std_s_safe = jnp.where(not_small, std_s, 1.0)

    C = jnp.where(not_small,
                  0.5 * (1.0 + jax.scipy.special.erf(s / (sqrt2 * std_s_safe))),
                  jnp.where(s > 0, 1.0, 0.0))
    qn = jnp.where(not_small & (C != 0.0),
                   s * C + (std_s / sqrt2pi)
                   * jnp.exp(-0.5 * (s / std_s_safe) ** 2),
                   jnp.where(~not_small & (s > 0), s, 0.0))

    qn_le_zero = qn <= 0.0
    C = jnp.where(qn_le_zero, 0.0, C)
    qn = jnp.where(qn_le_zero, 0.0, qn)
    return s, std_s, qn, C


def shoc_assumed_pdf_compute_sgs_liquid(a, ql1, ql2):
    """SGS liquid water: max(0, a ql1 + (1-a) ql2)."""
    return jnp.maximum(0.0, a * ql1 + (1.0 - a) * ql2)


def shoc_assumed_pdf_compute_cloud_liquid_variance(a, s1, ql1, C1, std_s1,
                                                   s2, ql2, C2, std_s2, shoc_ql):
    """Cloud liquid variance (CLUBB formulation adjusted to SHOC params)."""
    return jnp.maximum(0.0,
                       a * (s1 * ql1 + C1 * std_s1 ** 2)
                       + (1.0 - a) * (s2 * ql2 + C2 * std_s2 ** 2)
                       - shoc_ql ** 2)


def shoc_assumed_pdf_compute_liquid_water_flux(a, w1_1, w_first, ql1, w1_2, ql2):
    """Liquid water flux: a (w1_1-w) ql1 + (1-a)(w1_2-w) ql2."""
    return a * ((w1_1 - w_first) * ql1) + (1.0 - a) * ((w1_2 - w_first) * ql2)


def shoc_assumed_pdf_compute_buoyancy_flux(wthlsec, wqwsec, pval, wqls):
    """SGS buoyancy flux (Functions::shoc_assumed_pdf_compute_buoyancy_flux)."""
    epsterm = c.Rair / c.RV
    return (wthlsec + ((1.0 - epsterm) / epsterm) * c.basetemp * wqwsec
            + ((c.LatVap / c.CP) * (c.P0 / pval) ** (c.Rair / c.CP)
               - (1.0 / epsterm) * c.basetemp) * wqls)


@functools.partial(jax.jit, static_argnames=("extra_diags",))
def shoc_assumed_pdf(thetal, qw, w_field, thl_sec, qw_sec, dtime,
                     extra_diags: bool, wthl_sec, w_sec, wqw_sec, qwthl_sec,
                     w3, pres, zt_grid, zi_grid, shoc_ql_in):
    """Assumed-PDF driver (Functions::shoc_assumed_pdf).

    Midpoint inputs: thetal, qw, w_field, w_sec, pres, shoc_ql_in (previous
    ql, used for the cond/evap diagnostics). Interface inputs: thl_sec,
    qw_sec, wthl_sec, wqw_sec, qwthl_sec, w3.

    Returns (shoc_cldfrac, shoc_ql, wqls, wthv_sec, shoc_ql2,
    shoc_cond, shoc_evap); the last two are zeros unless extra_diags.
    """
    zt_grid = jnp.asarray(zt_grid)
    thetal = jnp.asarray(thetal)
    nlev = zt_grid.shape[-1]
    nlevi = nlev + 1
    zi_grid = jnp.asarray(zi_grid)

    # Interface -> midpoint interpolation (floors as in the C++)
    w3_zt = linear_interp(zi_grid, zt_grid, jnp.asarray(w3), nlevi, nlev, sc.largeneg)
    thl_sec_zt = linear_interp(zi_grid, zt_grid, jnp.asarray(thl_sec), nlevi, nlev, 0.0)
    wthl_sec_zt = linear_interp(zi_grid, zt_grid, jnp.asarray(wthl_sec), nlevi, nlev, sc.largeneg)
    qwthl_sec_zt = linear_interp(zi_grid, zt_grid, jnp.asarray(qwthl_sec), nlevi, nlev, sc.largeneg)
    wqw_sec_zt = linear_interp(zi_grid, zt_grid, jnp.asarray(wqw_sec), nlevi, nlev, sc.largeneg)
    qw_sec_zt = linear_interp(zi_grid, zt_grid, jnp.asarray(qw_sec), nlevi, nlev, 0.0)

    # Previous ql, with the top level zeroed (transcribed from the C++,
    # which zeroes shoc_ql[level 0] before the diagnostics).
    shoc_ql_prev = jnp.asarray(shoc_ql_in).at[..., 0].set(0.0)

    w_sec_m = jnp.asarray(w_sec)
    pval = jnp.asarray(pres)
    qw_m = jnp.asarray(qw)
    w_first = jnp.asarray(w_field)

    sqrtw2 = jnp.sqrt(w_sec_m)
    sqrtthl = jnp.maximum(_THL_TOL, jnp.sqrt(thl_sec_zt))
    sqrtqt = jnp.maximum(_RT_TOL, jnp.sqrt(qw_sec_zt))

    skew_w, w1_1, w1_2, w2_1, w2_2, a = shoc_assumed_pdf_vv_parameters(
        w_first, w_sec_m, w3_zt)

    thl1_1, thl1_2, thl2_1, thl2_2, sqrtthl2_1, sqrtthl2_2 = \
        shoc_assumed_pdf_thl_parameters(
            wthl_sec_zt, sqrtw2, sqrtthl, thl_sec_zt, thetal,
            w1_1, w1_2, skew_w, a)

    qw1_1, qw1_2, qw2_1, qw2_2, sqrtqw2_1, sqrtqw2_2 = \
        shoc_assumed_pdf_qw_parameters(
            wqw_sec_zt, sqrtw2, skew_w, sqrtqt, qw_sec_zt,
            w1_2, w1_1, qw_m, a)

    w1_1 = shoc_assumed_pdf_tilde_to_real(w_first, sqrtw2, w1_1)
    w1_2 = shoc_assumed_pdf_tilde_to_real(w_first, sqrtw2, w1_2)

    r_qwthl_1 = shoc_assumed_pdf_inplume_correlations(
        sqrtqw2_1, sqrtthl2_1, a, sqrtqw2_2, sqrtthl2_2, qwthl_sec_zt,
        qw1_1, qw_m, thl1_1, thetal, qw1_2, thl1_2)

    Tl1_1 = shoc_assumed_pdf_compute_temperature(thl1_1, pval)
    Tl1_2 = shoc_assumed_pdf_compute_temperature(thl1_2, pval)
    Tl1_1 = jnp.where(Tl1_1 <= _TL_MIN, _TL_MIN, Tl1_1)
    Tl1_2 = jnp.where(Tl1_2 <= _TL_MIN, _TL_MIN, Tl1_2)

    qs1, beta1, qs2, beta2 = shoc_assumed_pdf_compute_qs(Tl1_1, Tl1_2, pval)

    s1, std_s1, qn1, C1 = shoc_assumed_pdf_compute_s(
        qw1_1, qs1, beta1, pval, thl2_1, qw2_1, sqrtthl2_1, sqrtqw2_1, r_qwthl_1)
    s2c, std_s2c, qn2c, C2c = shoc_assumed_pdf_compute_s(
        qw1_2, qs2, beta2, pval, thl2_2, qw2_2, sqrtthl2_2, sqrtqw2_2, r_qwthl_1)

    equal = (qw1_1 == qw1_2) & (thl2_1 == thl2_2) & (qs1 == qs2)
    s2 = jnp.where(equal, s1, s2c)
    std_s2 = jnp.where(equal, std_s1, std_s2c)
    qn2 = jnp.where(equal, qn1, qn2c)
    C2 = jnp.where(equal, C1, C2c)

    ql1 = jnp.minimum(qn1, qw1_1)
    ql2 = jnp.minimum(qn2, qw1_2)

    shoc_cldfrac = jnp.minimum(1.0, a * C1 + (1.0 - a) * C2)

    if extra_diags:
        dum = jnp.maximum(0.0, a * ql1 + (1.0 - a) * ql2)
        shoc_cond = jnp.maximum(0.0, (dum - shoc_ql_prev) / dtime)
        shoc_evap = jnp.maximum(0.0, (shoc_ql_prev - dum) / dtime)
    else:
        shoc_cond = jnp.zeros_like(thetal)
        shoc_evap = jnp.zeros_like(thetal)

    shoc_ql = shoc_assumed_pdf_compute_sgs_liquid(a, ql1, ql2)
    shoc_ql2 = shoc_assumed_pdf_compute_cloud_liquid_variance(
        a, s1, ql1, C1, std_s1, s2, ql2, C2, std_s2, shoc_ql)
    wqls = shoc_assumed_pdf_compute_liquid_water_flux(
        a, w1_1, w_first, ql1, w1_2, ql2)
    wthv_sec = shoc_assumed_pdf_compute_buoyancy_flux(
        wthl_sec_zt, wqw_sec_zt, pval, wqls)

    return shoc_cldfrac, shoc_ql, wqls, wthv_sec, shoc_ql2, shoc_cond, shoc_evap
