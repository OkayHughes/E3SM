"""Port of eam/src/physics/cam/wv_sat_methods.F90 + wv_saturation.F90.

PORT_NOTES
----------
- Four SVP schemes (indices as in the Fortran): 0 OldGoffGratch,
  1 GoffGratch (default), 2 MurphyKoop, 3 Bolton. Bolton has no ice
  formula (`wv_sat_svp_ice` falls through to the water formula).
- tboil = 373.16 K is wv_saturation's own parameter (NOT physconst);
  ttrice = 20 K transition range; table spans tmin=127.16 to
  tmax=375.16 K in 1 K steps (plenest = ceil(tmax-tmin)+2 = 250).
- `qsat` (mixed-phase) reads the table via estblf's linear
  interpolation; `qsat_water`/`qsat_ice` evaluate the default scheme
  directly. All three clamp qs to 1 when p <= es and then clamp
  es = min(es, p).
- Derivative outputs follow calc_hltalt: above tmelt latvap loses
  2369*(t-tmelt); below tmelt a weight*latice is added with weight
  ramping over ttrice, plus the pcf polynomial tterm for d(es)/dT in
  the transition range (mixed-phase only).
- findsp is an 8-iteration Newton with three exits (converged -> 0,
  t < tmin bailout -> 4, non-convergence -> 2) plus a status-1
  feasibility pre-check and a final unconditional enthalpy check that
  overwrites status with 8. The port replicates the exact update and
  exit order with masked lax-style iteration so results match the
  Fortran trajectory point by point.

Everything is float64 and shape-polymorphic (elementwise over arrays).
"""

import jax.numpy as jnp
import numpy as np

from .constants import (CPAIR, EPSILO, H2OTRIP, LATICE, LATVAP, RH2O,
                        TMELT)

TBOIL = 373.16
TTRICE = 20.0
TMIN = 127.16
TMAX = 375.16
OMEPS = 1.0 - EPSILO

# polynomial for d(es)/dT in the water->ice transition (wv_saturation pcf)
PCF = np.array([5.04469588506e-01, -5.47288442819e+00,
                -3.67471858735e-01, -8.95963532403e-03,
                -7.78053686625e-05])

OLDGOFFGRATCH_IDX = 0
GOFFGRATCH_IDX = 1
MURPHYKOOP_IDX = 2
BOLTON_IDX = 3
DEFAULT_IDX = GOFFGRATCH_IDX


# ---------------------------------------------------------------------------
# SVP schemes (wv_sat_methods)
# ---------------------------------------------------------------------------
def goffgratch_svp_water(t):
    tmp1 = -7.90298 * (TBOIL / t - 1.0)
    tmp2 = 5.02808 * jnp.log10(TBOIL / t)
    tmp3 = 1.3816e-7 * (10.0 ** (11.344 * (1.0 - t / TBOIL)) - 1.0)
    tmp4 = 8.1328e-3 * (10.0 ** (-3.49149 * (TBOIL / t - 1.0)) - 1.0)
    tmp5 = jnp.log10(1013.246)
    return 10.0 ** (tmp1 + tmp2 - tmp3 + tmp4 + tmp5) * 100.0


def goffgratch_svp_ice(t):
    return 10.0 ** (-9.09718 * (H2OTRIP / t - 1.0)
                    - 3.56654 * jnp.log10(H2OTRIP / t)
                    + 0.876793 * (1.0 - t / H2OTRIP)
                    + jnp.log10(6.1071)) * 100.0


def murphykoop_svp_water(t):
    return jnp.exp(54.842763 - (6763.22 / t) - (4.210 * jnp.log(t))
                   + (0.000367 * t)
                   + (jnp.tanh(0.0415 * (t - 218.8))
                      * (53.878 - (1331.22 / t) - (9.44523 * jnp.log(t))
                         + 0.014025 * t)))


def murphykoop_svp_ice(t):
    return jnp.exp(9.550426 - (5723.265 / t) + (3.53068 * jnp.log(t))
                   - (0.00728332 * t))


def oldgoffgratch_svp_water(t):
    ps = 1013.246
    e1 = 11.344 * (1.0 - t / TBOIL)
    e2 = -3.49149 * (TBOIL / t - 1.0)
    f1 = -7.90298 * (TBOIL / t - 1.0)
    f2 = 5.02808 * jnp.log10(TBOIL / t)
    f3 = -1.3816 * (10.0 ** e1 - 1.0) / 10000000.0
    f4 = 8.1328 * (10.0 ** e2 - 1.0) / 1000.0
    f5 = jnp.log10(ps)
    return (10.0 ** (f1 + f2 + f3 + f4 + f5)) * 100.0


def oldgoffgratch_svp_ice(t):
    term1 = 2.01889049 / (TMELT / t)
    term2 = 3.56654 * jnp.log(TMELT / t)
    term3 = 20.947031 * (TMELT / t)
    return 575.185606e10 * jnp.exp(-(term1 + term2 + term3))


def bolton_svp_water(t):
    c1, c2, c3 = 611.2, 17.67, 243.5
    return c1 * jnp.exp((c2 * (t - TMELT)) / ((t - TMELT) + c3))


_WATER = {OLDGOFFGRATCH_IDX: oldgoffgratch_svp_water,
          GOFFGRATCH_IDX: goffgratch_svp_water,
          MURPHYKOOP_IDX: murphykoop_svp_water,
          BOLTON_IDX: bolton_svp_water}
_ICE = {OLDGOFFGRATCH_IDX: oldgoffgratch_svp_ice,
        GOFFGRATCH_IDX: goffgratch_svp_ice,
        MURPHYKOOP_IDX: murphykoop_svp_ice,
        BOLTON_IDX: bolton_svp_water}  # Bolton has no ice formula


def svp_water(t, idx=DEFAULT_IDX):
    return _WATER[idx](jnp.asarray(t, dtype=jnp.float64))


def svp_ice(t, idx=DEFAULT_IDX):
    return _ICE[idx](jnp.asarray(t, dtype=jnp.float64))


def svp_trans(t, idx=DEFAULT_IDX):
    """Water/ice transition-weighted SVP (wv_sat_svp_trans)."""
    t = jnp.asarray(t, dtype=jnp.float64)
    es = jnp.where(t >= TMELT - TTRICE, svp_water(t, idx), 0.0)
    esice = svp_ice(t, idx)
    weight = jnp.where(TMELT - t > TTRICE, 1.0, (TMELT - t) / TTRICE)
    return jnp.where(t < TMELT, weight * esice + (1.0 - weight) * es, es)


def svp_to_qsat(es, p):
    """Saturation specific humidity from SVP; qs = 1 where p <= es."""
    return jnp.where(p - es <= 0.0, 1.0,
                     EPSILO * es / (p - OMEPS * es))


# ---------------------------------------------------------------------------
# Table (wv_saturation estbl / estblf)
# ---------------------------------------------------------------------------
def build_estbl():
    """SVP lookup table: svp_trans at tmin + i, i = 0..plenest-1."""
    plenest = int(np.ceil(TMAX - TMIN)) + 2
    return svp_trans(TMIN + jnp.arange(plenest, dtype=jnp.float64))


_ESTBL = None


def _estbl():
    global _ESTBL
    if _ESTBL is None:
        _ESTBL = build_estbl()
    return _ESTBL


def estblf(t, estbl=None):
    """Linear interpolation in the SVP table (Fortran 1-based indexing
    becomes 0-based here)."""
    if estbl is None:
        estbl = _estbl()
    t = jnp.asarray(t, dtype=jnp.float64)
    t_tmp = jnp.maximum(jnp.minimum(t, TMAX) - TMIN, 0.0)
    i = t_tmp.astype(jnp.int64)          # int() truncation, as in Fortran
    weight = t_tmp - jnp.trunc(t_tmp)    # aint()
    return (1.0 - weight) * estbl[i] + weight * estbl[i + 1]


# ---------------------------------------------------------------------------
# Latent-heat corrections and derivative outputs
# ---------------------------------------------------------------------------
def no_ip_hltalt(t):
    hltalt = jnp.full_like(t, LATVAP)
    return jnp.where(t >= TMELT, hltalt - 2369.0 * (t - TMELT), hltalt)


def calc_hltalt(t):
    """Returns (hltalt, tterm) — ice-phase-aware latent heat and the
    transition-region d(es)/dT polynomial term."""
    hltalt = no_ip_hltalt(t)
    tc = t - TMELT
    in_trans = tc >= -TTRICE
    weight = jnp.where(in_trans, -tc / TTRICE, 1.0)
    # Horner evaluation of pcf, exactly as the Fortran loop
    tterm = jnp.zeros_like(t)
    for i in range(len(PCF) - 1, -1, -1):
        tterm = PCF[i] + tc * tterm
    tterm = jnp.where(in_trans & (t < TMELT), tterm / TTRICE, 0.0)
    hltalt = jnp.where(t < TMELT, hltalt + weight * LATICE, hltalt)
    return hltalt, tterm


def tq_enthalpy(t, q, hltalt):
    return CPAIR * t + hltalt * q


def deriv_outputs(t, p, es, qs, hltalt, tterm):
    """(gam, dqsdt) as in wv_saturation deriv_outputs."""
    desdt = hltalt * es / (RH2O * t * t) + tterm
    dqsdt = jnp.where(qs == 1.0, 0.0,
                      qs * p * desdt / (es * (p - OMEPS * es)))
    return dqsdt * (hltalt / CPAIR), dqsdt


# ---------------------------------------------------------------------------
# qsat family — return dict(es, qs, gam, dqsdt, enthalpy)
# ---------------------------------------------------------------------------
def _finish(t, p, es, qs, hltalt, tterm):
    gam, dqsdt = deriv_outputs(t, p, es, qs, hltalt, tterm)
    return {"es": es, "qs": qs, "gam": gam, "dqsdt": dqsdt,
            "enthalpy": tq_enthalpy(t, qs, hltalt)}


def qsat(t, p, estbl=None):
    """Mixed-phase (table) saturation, with derivatives."""
    t = jnp.asarray(t, dtype=jnp.float64)
    p = jnp.asarray(p, dtype=jnp.float64)
    es = estblf(t, estbl)
    qs = svp_to_qsat(es, p)
    es = jnp.minimum(es, p)
    hltalt, tterm = calc_hltalt(t)
    return _finish(t, p, es, qs, hltalt, tterm)


def qsat_water(t, p, idx=DEFAULT_IDX):
    t = jnp.asarray(t, dtype=jnp.float64)
    p = jnp.asarray(p, dtype=jnp.float64)
    es = svp_water(t, idx)
    qs = svp_to_qsat(es, p)
    es = jnp.minimum(es, p)
    hltalt = no_ip_hltalt(t)
    return _finish(t, p, es, qs, hltalt, jnp.zeros_like(t))


def qsat_ice(t, p, idx=DEFAULT_IDX):
    t = jnp.asarray(t, dtype=jnp.float64)
    p = jnp.asarray(p, dtype=jnp.float64)
    es = svp_ice(t, idx)
    qs = svp_to_qsat(es, p)
    es = jnp.minimum(es, p)
    hltalt = jnp.full_like(t, LATVAP + LATICE)
    return _finish(t, p, es, qs, hltalt, jnp.zeros_like(t))


# ---------------------------------------------------------------------------
# findsp — wet-bulb temperature Newton solver
# ---------------------------------------------------------------------------
C3 = 287.04 * (7.5 * np.log(10.0)) / CPAIR  # wv_sat_init's hardcoded 287.04


def _qsat_choice(t, p, use_ice, estbl):
    return qsat(t, p, estbl) if use_ice else qsat_water(t, p)


def _hltalt_choice(t, use_ice):
    return calc_hltalt(t)[0] if use_ice else no_ip_hltalt(t)


def findsp(q, t, p, use_ice, estbl=None):
    """Wet-bulb temperature/humidity (elemental findsp), vectorized.

    Returns (tsp, qsp, status) with the Fortran status codes. The
    8-iteration Newton loop is replicated with masked updates so each
    point follows exactly the Fortran trajectory (including freezing
    its state on convergence or the tmin bailout).
    """
    q = jnp.asarray(q, dtype=jnp.float64)
    t = jnp.asarray(t, dtype=jnp.float64)
    p = jnp.asarray(p, dtype=jnp.float64)
    if estbl is None and use_ice:
        estbl = _estbl()

    r0 = _qsat_choice(t, p, use_ice, estbl)
    es0, qs0 = r0["es"], r0["qs"]

    unusable = ((p <= 5.0 * es0) | (qs0 <= 0.0) | (qs0 >= 0.5)
                | (t < TMIN) | (t > TMAX))

    hltalt0 = _hltalt_choice(t, use_ice)
    enin = tq_enthalpy(t, q, hltalt0)

    # UKMO first guess
    c1 = hltalt0 * C3
    c2 = (t + 36.0) ** 2
    r1b = c2 / (c2 + c1 * qs0)
    qvd = r1b * (q - qs0)
    tsp = t + (hltalt0 / CPAIR) * qvd

    rg = _qsat_choice(tsp, p, use_ice, estbl)
    qsp, gam, enout = rg["qs"], rg["gam"], rg["enthalpy"]

    status = jnp.full(t.shape, 2, dtype=jnp.int64)
    done = jnp.zeros(t.shape, dtype=bool)

    dttol = 1.0e-4
    dqtol = 1.0e-4

    for _ in range(8):
        active = ~done
        g = enin - enout
        dgdt = -CPAIR * (1.0 + gam)
        t1 = tsp - g / dgdt
        dt_rel = jnp.abs(t1 - tsp) / t1
        tsp_new = t1

        # tmin bailout branch
        bail = tsp_new < TMIN
        tsp_bail = jnp.full_like(tsp_new, TMIN)
        hlt_bail = _hltalt_choice(tsp_bail, use_ice)
        qsp_bail = (enin - CPAIR * tsp_bail) / hlt_bail
        enout_bail = tq_enthalpy(tsp_bail, qsp_bail, hlt_bail)

        # normal branch: re-evaluate at tsp_new
        rn = _qsat_choice(tsp_new, p, use_ice, estbl)
        q1, gam_new, enout_new = rn["qs"], rn["gam"], rn["enthalpy"]
        dq_rel = jnp.abs(q1 - qsp) / jnp.maximum(q1, 1.0e-12)
        conv = (dt_rel < dttol) & (dq_rel < dqtol)

        upd = lambda new, old: jnp.where(active, new, old)  # noqa: E731
        tsp = upd(jnp.where(bail, tsp_bail, tsp_new), tsp)
        qsp = upd(jnp.where(bail, qsp_bail, q1), qsp)
        gam = upd(jnp.where(bail, gam, gam_new), gam)
        enout = upd(jnp.where(bail, enout_bail, enout_new), enout)
        status = upd(jnp.where(bail, 4,
                               jnp.where(conv, 0, status)), status)
        done = done | (active & (bail | conv))

    # final enthalpy check overwrites status
    bad_enthalpy = jnp.abs((enin - enout) / (enin + enout)) > 1.0e-4
    status = jnp.where(bad_enthalpy, 8, status)

    # status-1 points keep their inputs
    tsp = jnp.where(unusable, t, tsp)
    qsp = jnp.where(unusable, q, qsp)
    status = jnp.where(unusable, 1, status)
    return tsp, qsp, status
