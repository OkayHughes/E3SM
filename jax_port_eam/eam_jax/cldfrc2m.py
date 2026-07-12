"""Port of the EAMv3 RH-based stratiform cloud-fraction module,
eam/src/physics/cam/cldfrc2m.F90: astG_PDF (triangular-PDF liquid
stratus fraction a and inverse gradient G = dU/da), astG_RHU (CAM3.5
quadratic formula) and aist (ice stratus fraction, iceopt 1-7).

PORT_NOTES
----------
- Everything is elementwise over arrays of any shape; the Fortran
  `_single` and vector routines contain byte-identical per-column
  bodies, so each pair is ported once (goldens replay both driver
  variants against the same function). One asymmetry: aist_single has
  NO iceopt=6 branch (its `aist` would be left undefined); this port
  implements the vector semantics, which are identical for every other
  iceopt. EAMv3 macrophysics only ever uses iceopt=5.
- Parameters: cldfrc2m_init copies rhminl/rhminl_adj_land/rhminh/
  premit/premib/iceopt/icecrit/minice out of cloud_fraction's
  cldfrc_readnl; rhmini/rhmaxi come from the cldfrc2m_nl namelist.
  make_params defaults are the EAMv3 phys="default" values from
  bld/namelist_files/namelist_defaults_eam.xml: rhminl=0.950
  (microphys="p3"), rhminl_adj_land=0.100, rhminh=0.800,
  premit=25000 Pa (dyn="se"), premib=70000 Pa, iceopt=5, icecrit=0.93,
  minice=1.0e-12, rhmini=0.80, rhmaxi=1.05 (clubb_sgs="1",
  clubb_do_deep="0").
- Compile-time module parameters, verbatim: qist_min=1e-7,
  qist_max=5e-3 (in-stratus IWC bounds), CAMstfrac=.false. (callers
  use astG_PDF, not astG_RHU), freeze_dry=.false. (the freeze-dry qv
  factor in the p>=premib band is dead code; ported behind the same
  constant), cldrh=1.0, mincld=1e-4.
- Pressure bands (all three schemes share them): p >= premib uses
  rhminl (minus rhminl_adj_land over snow-free land, snowh <= 1e-6 and
  nint(landfrac)==1); p < premit uses rhminh; between, rhmin blends
  with weight (premib - max(p,premit))/(premib - premit) and — as in
  the Fortran, where the land adjustment is commented out — ignores
  land.
- nint(landfrac)==1 rounds half AWAY from zero (landfrac=0.5 is land);
  jnp.round rounds half to even, so land is floor(landfrac + 0.5)==1
  (landfrac is in [0,1]).
- astG_PDF keeps the Fortran's literal 3.141592 (not pi) and the
  Ga=1e10 sentinels at a=0/a=1; astG_RHU's sentinel is 1e20 in the low
  and high bands but 1e10 in the mid band (verbatim quirk). Ga can be
  +inf inside the mid branch expressions when U==rhmin, but those
  lanes are masked to the sentinel exactly as the Fortran if-guard
  does.
- aist iceopt branches: 1 Wang&Sassen, 2 Schiller, 3 Wood&Field,
  4 Wilson&Ballard (literal 3.1415927), 5 modified Slingo (EAMv3
  default; RH over ice from (qv+qi)/qs * esl/esi with qsat_water /
  svp_water / svp_ice from eam_jax.wv_sat, GoffGratch default scheme),
  6 Gettelman/Heymsfield (vector-only), 7 all-or-nothing. For iceopt
  5/6 the minice/mincld limiter runs, then the icimr in-stratus bounds
  (icimr computed ONCE from the mincld-floored aist and tested against
  both qist_min and qist_max). Final clamp to [0, 0.999] everywhere.
- Real-exponent powers (** 2._r8, ** (2._r8/3._r8), bs**ttmp, nil**ch)
  are kept as float-exponent jnp powers; integer squares (**2) are
  x*x. Fortran/XLA libm last-ulp differences are inside the 1e-12
  replay tolerance.
"""

import jax.numpy as jnp

from .constants import RAIR
from .wv_sat import qsat_water, svp_ice, svp_water

# ---- module parameters (verbatim from cldfrc2m.F90) ----
QIST_MIN = 1.0e-7       # minimum in-stratus IWC [kg/kg]
QIST_MAX = 5.0e-3       # maximum in-stratus IWC [kg/kg]
CAMSTFRAC = False       # use astG_RHU instead of astG_PDF
FREEZE_DRY = False      # Vavrus & Waliser freeze-dry factor
CLDRH = 1.0             # in-stratus RH
MINCLD = 1.0e-4         # minimum ice cloud fraction (iceopt 5/6)
SNOWH_THRESH = 0.000001  # snow-free threshold [m lwe]

# aist fit constants (verbatim)
A_WS, B_WS, C_WS = 26.87, 0.569, 0.002892       # Wang & Sassen (1)
AS_SCH, BS_SCH, CS_SCH = -68.4202, 0.983917, 2.81795  # Schiller (2)
KC = 75.0                                       # Wood & Field (3)
AH, BH, CH = 6.73834e-08, 0.0533110, 0.3493813  # Gettelman (6)
PI_PDF = 3.141592       # astG_PDF's literal pi
PI_AIST = 3.1415927     # iceopt 4's literal pi


def make_params(rhminl=0.950, rhminl_adj_land=0.100, rhminh=0.800,
                premit=25000.0, premib=70000.0, iceopt=5,
                icecrit=0.93, minice=1.0e-12, rhmini=0.80,
                rhmaxi=1.05):
    """Runtime parameters; defaults are EAMv3 phys="default" (see
    PORT_NOTES)."""
    return dict(rhminl=rhminl, rhminl_adj_land=rhminl_adj_land,
                rhminh=rhminh, premit=premit, premib=premib,
                iceopt=int(iceopt), icecrit=icecrit, minice=minice,
                rhmini=rhmini, rhmaxi=rhmaxi)


def _land(landfrac):
    """Fortran nint(landfrac) == 1 (half away from zero; landfrac in
    [0, 1])."""
    return jnp.floor(jnp.asarray(landfrac, dtype=jnp.float64)
                     + 0.5) == 1.0


def _rhmin_liquid(p, landfrac, snowh, prm):
    """Critical RH of the liquid-stratus schemes: three pressure bands
    with the snow-free-land adjustment only in the lowest."""
    low = p >= prm["premib"]
    high = p < prm["premit"]
    land_sf = _land(landfrac) & (snowh <= SNOWH_THRESH)
    rhmin_low = jnp.where(land_sf,
                          prm["rhminl"] - prm["rhminl_adj_land"],
                          prm["rhminl"])
    rhwght = ((prm["premib"] - jnp.maximum(p, prm["premit"]))
              / (prm["premib"] - prm["premit"]))
    rhmin_mid = (prm["rhminh"] * rhwght
                 + prm["rhminl"] * (1.0 - rhwght))
    return jnp.where(low, rhmin_low,
                     jnp.where(high, prm["rhminh"], rhmin_mid))


def astG_PDF(U, p, qv, landfrac, snowh, prm):
    """Triangular-PDF stratus fraction. Returns (a, Ga, rhmin);
    rhmin is the Fortran orhmin output."""
    U = jnp.asarray(U, dtype=jnp.float64)
    p = jnp.asarray(p, dtype=jnp.float64)
    qv = jnp.asarray(qv, dtype=jnp.float64)
    snowh = jnp.asarray(snowh, dtype=jnp.float64)

    rhmin = _rhmin_liquid(p, landfrac, snowh, prm)
    dV = CLDRH - rhmin

    sqrt2 = jnp.sqrt(2.0)
    # branch 2: cldrh - dV/6 < U < 1
    a2 = 1.0 - (-3.0 / sqrt2 * (U - CLDRH) / dV) ** (2.0 / 3.0)
    ga2 = dV / sqrt2 * jnp.sqrt(1.0 - a2)
    # branch 3: cldrh - dV < U <= cldrh - dV/6
    a3 = 4.0 * jnp.cos((1.0 / 3.0)
                       * (jnp.arccos((3.0 / 2.0 / sqrt2)
                                     * (1.0 + (U - CLDRH) / dV))
                          - 2.0 * PI_PDF)) ** 2.0
    ga3 = dV / sqrt2 * (1.0 / jnp.sqrt(a3) - jnp.sqrt(a3))

    b1 = U >= 1.0
    b2 = (U > CLDRH - dV / 6.0) & (U < 1.0)
    b3 = (U > CLDRH - dV) & (U <= CLDRH - dV / 6.0)
    a = jnp.where(b1, 1.0,
                  jnp.where(b2, a2, jnp.where(b3, a3, 0.0)))
    ga = jnp.where(b1, 1.0e10,
                   jnp.where(b2, ga2, jnp.where(b3, ga3, 1.0e10)))

    if FREEZE_DRY:  # compile-time .false. in EAM (kept verbatim)
        fac = jnp.maximum(0.15, jnp.minimum(1.0, qv / 0.0030))
        low = p >= prm["premib"]
        a = jnp.where(low, a * fac, a)
        ga = jnp.where(low, ga / fac, ga)
    return a, ga, rhmin


def astG_RHU(U, p, qv, landfrac, snowh, prm):
    """CAM3.5 quadratic stratus fraction. Returns (a, Ga, rhmin)."""
    U = jnp.asarray(U, dtype=jnp.float64)
    p = jnp.asarray(p, dtype=jnp.float64)
    qv = jnp.asarray(qv, dtype=jnp.float64)
    snowh = jnp.asarray(snowh, dtype=jnp.float64)

    rhmin = _rhmin_liquid(p, landfrac, snowh, prm)
    rhdif = (U - rhmin) / (1.0 - rhmin)
    m = jnp.maximum(rhdif, 0.0)
    a = jnp.minimum(1.0, m * m)

    # sentinel is 1e20 in the low/high bands but 1e10 in the mid band
    low = p >= prm["premib"]
    high = p < prm["premit"]
    sentinel = jnp.where(low | high, 1.0e20, 1.0e10)
    ga_mid = 0.5 * (1.0 - rhmin) * ((1.0 - rhmin) / (U - rhmin))
    ga = jnp.where((U >= 1.0) | (U <= rhmin), sentinel, ga_mid)

    if FREEZE_DRY:  # compile-time .false. in EAM (kept verbatim)
        fac = jnp.maximum(0.15, jnp.minimum(1.0, qv / 0.0030))
        a = jnp.where(low, a * fac, a)
        ga = jnp.where(low, ga / fac, ga)
    return a, ga, rhmin


def aist(qv, T, p, qi, ni, landfrac, snowh, prm):
    """Ice stratus fraction (aist_vector semantics; see PORT_NOTES on
    aist_single/iceopt 6). ni is only used by iceopt=6."""
    qv = jnp.asarray(qv, dtype=jnp.float64)
    T = jnp.asarray(T, dtype=jnp.float64)
    p = jnp.asarray(p, dtype=jnp.float64)
    qi = jnp.asarray(qi, dtype=jnp.float64)
    ni = jnp.asarray(ni, dtype=jnp.float64)
    snowh = jnp.asarray(snowh, dtype=jnp.float64)
    iceopt = prm["iceopt"]

    qs = qsat_water(T, p)["qs"]
    esl = svp_water(T)
    esi = svp_ice(T)

    if iceopt < 3:
        if iceopt == 1:
            ttmp = jnp.maximum(195.0, jnp.minimum(T, 253.0)) - 273.16
            icicval = (A_WS + B_WS * ttmp
                       + C_WS * jnp.power(ttmp, 2.0))
            rho = p / (RAIR * T)
            icicval = icicval * 1.0e-6 / rho
        else:
            ttmp = jnp.maximum(190.0, jnp.minimum(T, 273.16))
            icicval = 10.0 ** (AS_SCH * jnp.power(BS_SCH, ttmp)
                               + CS_SCH)
            icicval = icicval * 1.0e-6 * 18.0 / 28.97
        ist = jnp.maximum(0.0, jnp.minimum(qi / icicval, 1.0))
    elif iceopt == 3:
        ist = 1.0 - jnp.exp(-KC * qi / (qs * (esi / esl)))
        ist = jnp.maximum(0.0, jnp.minimum(ist, 1.0))
    elif iceopt == 4:
        ncf = qi / ((1.0 - prm["icecrit"]) * qs)
        phi = (jnp.arccos(jnp.clip(3.0 * (1.0 - ncf) / 2.0 ** 1.5,
                                   -1.0, 1.0))
               + 4.0 * PI_AIST) / 3.0
        # clip only protects the masked-out lanes (|arg|>1 implies the
        # branch below is not selected); in-branch values are untouched
        ist_mid = 1.0 - 4.0 * jnp.cos(phi) * jnp.cos(phi)
        ist_lo = 0.5 * jnp.power(6.0 * ncf, 2.0 / 3.0)
        ist = jnp.where(
            ncf <= 0.0, 0.0,
            jnp.where(ncf <= 1.0 / 6.0, ist_lo,
                      jnp.where(ncf < 1.0, ist_mid, 1.0)))
        ist = jnp.maximum(0.0, jnp.minimum(ist, 1.0))
    elif iceopt == 5:
        rhi = (qv + qi) / qs * (esl / esi)
        rhdif = (rhi - prm["rhmini"]) / (prm["rhmaxi"] - prm["rhmini"])
        m = jnp.maximum(rhdif, 0.0)
        ist = jnp.minimum(1.0, m * m)
    elif iceopt == 6:
        rho = p / (RAIR * T)
        nil = ni * rho / 1000.0
        icicval = AH * jnp.exp(BH * T) * jnp.power(nil, CH)
        icicval = icicval / rho / 1000.0
        ist = jnp.maximum(0.0, jnp.minimum(qi / icicval, 1.0))
        ist = jnp.minimum(ist, 1.0)
    elif iceopt == 7:
        ist = jnp.where(qi >= prm["minice"], 1.0, 0.0)
    else:
        raise ValueError(f"unsupported iceopt {iceopt}")

    if iceopt in (5, 6):
        # empty-cloud / no-cloud-ice limiter, then in-stratus IWC
        # bounds; icimr evaluated once from the mincld-floored value
        ist = jnp.where(qi < prm["minice"], 0.0,
                        jnp.maximum(MINCLD, ist))
        has_ice = qi >= prm["minice"]
        icimr = qi / jnp.where(ist > 0.0, ist, 1.0)  # guarded: ist >=
        # mincld wherever has_ice, and icimr is only used under has_ice
        ist = jnp.where(has_ice & (icimr < QIST_MIN),
                        jnp.maximum(0.0,
                                    jnp.minimum(1.0, qi / QIST_MIN)),
                        ist)
        ist = jnp.where(has_ice & (icimr > QIST_MAX),
                        jnp.maximum(0.0,
                                    jnp.minimum(1.0, qi / QIST_MAX)),
                        ist)

    # 0.999 cap prevents infinite ql_st in instratus_condensate
    return jnp.maximum(0.0, jnp.minimum(ist, 0.999))
