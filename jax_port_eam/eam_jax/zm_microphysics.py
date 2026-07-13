"""Port of the ZM convective microphysics,
eam/src/physics/cam/zm/zm_microphysics.F90 (zm_mphy + zm_mphyi +
zm_microphysics_adjust), together with its aerosol-activation
dependencies cam/activate_drop_mam.F90 (actdrop_mam_init/_calc,
Abdul-Razzak & Ghan 2000) and cam/nucleate_ice_conv.F90
(nucleati_conv, Liu & Penner 2005).

PORT_NOTES
----------
- ACTIVATION STRATEGY (PORTING_PLAN.md row 9): the real activation
  modules compile cleanly against the existing infrastructure stubs,
  so the Tier-1 harness (harness/build_zm_microp.py) compiles them
  UNMODIFIED and this module ports them outright - no prescribed
  activation values anywhere. Only the BULK aerosol scheme
  (ndrop_bam) is out of scope: EAMv3 production runs MAM modal
  aerosols, the harness stubs ndrop_bam abort-only, and this port
  raises on aero["scheme"] != "modal".
- actdrop_mam_calc is ported for the exact configuration zm_mphy
  uses: sigw = wmix = 0 -> the SINGLE-UPDRAFT branch only (the
  spectral-updraft integration loop is never reached; documented
  scope). Of its outputs only fn feeds back (nlsrc); fm/fluxn/fluxm/
  flux_fullact are not computed. maxsat's mode scan ("all modes weak
  -> smax=1e-20, else sum over all modes") is reproduced exactly.
- CPP variant: compiled without MODAL_AERO_4MODE_MOM /
  MODAL_AERO_5MODE / RAIN_EVAP_TO_COARSE_AERO, so the coarse-mode
  dust weight is wght = dmc/(ssmc + dmc + so4mc) and the aero object
  carries only dust/nacl/so4 coarse species indices.
- gamma: with gfortran, shr_spfn_gamma resolves to the compiler's
  gamma intrinsic (glibc tgamma) via HAVE_GAMMA_INTRINSICS; the port
  uses scipy.special.gamma (cephes) on concrete values - the two
  libms agree to ~1 ulp, which is inside the measured replay
  tolerances. Same story for erf (glibc erf vs XLA erf). Constant
  gamma arguments (zm_mphyi cons*, fallspeed factors) are evaluated
  once at init in float64.
- Single-precision literal promotion reproduced verbatim: mg0 =
  1.6E-10 (default-real parameter), the qric/qgic >= 1.e-8
  thresholds in the rain-graupel collection branch, and the
  (rhosu/rho)**0.54 exponents in the two POST-integration fallspeed
  blocks (all other 0.54 exponents are _r8) are float(np.float32(x)).
- Index convention: 0-based, level 0 = model top, exactly as
  eam_jax/zm_conv.py: Fortran midpoint k maps to python k-1; msg
  keeps its Fortran VALUE; jb/jt/jlcl inputs are 0-based. The Fortran
  main loop do k = pver, msg+2, -1 is python kk = pver-1 .. msg+1.
- Structure: the Fortran column/iteration/level loops (i, it=1..2,
  k descending) are reordered to it-outer, k-outer, columns
  vectorized with per-column masks - exact, since columns are
  independent. Every sequential dependence (vertical integration
  writing level k-1, the it=2 "falling from above" sums that read
  it=1 fallspeeds/states of levels above, kqc/kqi boundary flags) is
  preserved by keeping full (ncol, pver) arrays and replaying the
  exact statement order. Reductions accumulate in Fortran loop order.
- Faithful quirks kept:
  * columns glaciated from cloud base (cmei > qsmall at the first
    plume level) never trigger the kqi ice boundary condition, so
    qi/ni stay 0 while frz still carries cmei (goldened, config k
    "ice" family);
  * npra = pra/(qcic/ncic) etc. rely on IEEE (x/0=inf, y/inf=0);
    reproduced with the same expressions on masked-safe lanes;
  * the bulk-scheme immersion-freezing branch can read stale
    mnuccc/nnuccc elements - unreachable in the modal scope;
  * nsagg is "hard-wired for bs=0.4" via generic pow expressions -
    transcribed verbatim;
  * uninitialized-diagnostic subtleties (nihf/niimm/nidep/nimey are
    write-only locals) are dropped: they never feed back.
- lamc/pgam are intent(inout): seeded by the caller with
  loc_microp_st%lambdadpcu = (mucon+1)/dcon and mudpcu = mucon.
- deltat comes from get_step_size() in the Fortran (time_manager
  stub); here it is an explicit argument, as are auto_fac/accr_fac/
  dcs (zmconv_auto_fac/accr_fac/micro_dcs) and grav/cp/rd (zm_const).
- Everything float64, eager-mode (concrete masks and python loops;
  not jit-compatible as written, matching eam_jax/zm_conv.py).
"""

import numpy as np
import jax.numpy as jnp
from scipy.special import gamma as _scipy_gamma
from jax.scipy.special import erf as _erf

from .constants import (GRAVIT, RAIR, TMELT, CPAIR, RH2O, RGAS, MWH2O,
                        RHOH2O, LATVAP, LATICE, EPSILO, SHR_CONST_PI)
from .wv_sat import svp_water, svp_ice, qsat as _qsat_pa

__all__ = ["make_mphyi", "make_actdrop_params", "actdrop_mam_calc",
           "nucleati_conv", "zm_mphy", "zm_microphysics_adjust",
           "MUCON", "DCON", "QSMALL"]

W = jnp.where

# ---------------------------------------------------------------------------
# module constants (zm_microphysics.F90 header, verbatim)
# ---------------------------------------------------------------------------
PI = 3.14159265358979323846
RHOW, RHOI, RHOG, RHOSN = 1000.0, 500.0, 400.0, 100.0
AC, BC = 3.0e7, 2.0        # droplet fallspeed V = a D^b
AS_, BS = 11.72, 0.41      # snow
AI, BI = 700.0, 1.0        # cloud ice
AR, BR = 841.99667, 0.8    # rain
AG, BG = 19.3, 0.37        # graupel
EII = 0.1
ECR = 1.0
ECG = 0.7
BIMM, AIMM = 100.0, 0.66
MG0 = float(np.float32(1.6e-10))   # default-real literal in the Fortran
RIN = 0.1e-6
KA_B, LS_B, RV_B = 2.4e-2, 2.834e6, 461.0
RN_DST1, RN_DST2 = 0.258e-6, 0.717e-6
RN_DST3, RN_DST4 = 1.576e-6, 3.026e-6
QSMALL = 1.0e-18
DCON = 25.0e-6
MUCON = 5.3
LAMBDADPCU = (MUCON + 1.0) / DCON
# single-precision literals promoted to double (see PORT_NOTES)
_THR_1EM8_SP = float(np.float32(1.0e-8))
_EXP_054_SP = float(np.float32(0.54))


def _gamma(x):
    """shr_spfn_gamma: with gfortran it is the gamma intrinsic (glibc
    tgamma); scipy's cephes gamma matches to ~1 ulp. Eager-only."""
    return jnp.asarray(_scipy_gamma(np.asarray(x, dtype=np.float64)))


def _divs(a, s):
    """Fortran-faithful division of an array by a scalar. XLA rewrites
    x / const into x * (1/const), which is NOT correctly rounded and
    flips razor-edge branch comparisons; materializing the divisor
    forces a true IEEE division. (Divisions by powers of two are exact
    either way and stay plain.)"""
    return a / jnp.full_like(a, s)


def _bdiv(a, b):
    """Fortran-faithful division by a broadcast array divisor (XLA
    applies the same reciprocal rewrite to broadcast divisors)."""
    return a / jnp.broadcast_to(b, a.shape)


def make_mphyi():
    """zm_mphyi: derived microphysics constants."""
    g = lambda x: float(_scipy_gamma(np.float64(x)))  # noqa: E731
    d = {}
    d["xlf"] = LATICE
    d["ci"], d["di"] = RHOI * PI / 6.0, 3.0
    d["cs"], d["ds"] = RHOSN * PI / 6.0, 3.0
    d["cr"], d["dr"] = RHOW * PI / 6.0, 3.0
    d["cg"], d["dg"] = RHOG * PI / 6.0, 3.0
    d["rhosu"] = 85000.0 / (RAIR * TMELT)
    d["alfa_b"] = 1.0 / 3.0
    d["rhoi13"] = RHOI ** d["alfa_b"]
    d["c23"] = 2.0 / 3.0
    d["mi0"] = 4.0 / 3.0 * PI * RHOI * (10.0e-6) * (10.0e-6) * (10.0e-6)
    d["mmult"] = 4.0 / 3.0 * PI * RHOI * (5.0e-6) ** 3
    d["cons14"] = g(BG + 3.0) * PI / 4.0 * ECG
    d["cons16"] = g(BI + 3.0) * PI / 4.0 * ECG
    d["cons17"] = (4.0 * 2.0 * 3.0 * d["rhosu"] * PI * ECG * ECG
                   * g(2.0 * BS + 2.0) / (8.0 * (RHOG - RHOSN)))
    d["cons18"] = RHOSN * RHOSN
    d["cons19"] = RHOW * RHOW
    d["cons24"] = PI / 4.0 * ECR * g(BR + 3.0)
    d["cons25"] = PI * PI / 24.0 * RHOW * ECR * g(BR + 6.0)
    d["cons31"] = PI * PI * ECR * RHOSN
    d["cons32"] = PI / 2.0 * ECR
    d["cons41"] = PI * PI * ECR * RHOW
    d["droplet_mass_25um"] = 4.0 / 3.0 * PI * RHOW * (25.0e-6) ** 3
    # constant gamma factors used repeatedly in zm_mphy
    d["g_1di"] = g(1.0 + d["di"])
    d["g_1ds"] = g(1.0 + d["ds"])
    d["g_1dg"] = g(1.0 + d["dg"])
    d["g_1br"] = g(1.0 + BR)
    d["g_4br"] = g(4.0 + BR)
    d["g_1bs"] = g(1.0 + BS)
    d["g_4bs"] = g(4.0 + BS)
    d["g_1bg"] = g(1.0 + BG)
    d["g_4bg"] = g(4.0 + BG)
    d["g_bs3"] = g(BS + 3.0)
    return d


# ---------------------------------------------------------------------------
# activate_drop_mam (Abdul-Razzak & Ghan 2000), single-updraft branch
# ---------------------------------------------------------------------------
T0_ACT = 273.0
P0_ACT = 1013.25e2
SURFTEN = 0.076
THIRD = 1.0 / 3.0
TWOTHIRD = 2.0 * THIRD
SQ2 = float(np.sqrt(np.float64(2.0)))


def make_actdrop_params(sigmag_amode):
    """actdrop_mam_init: per-mode width functions + aten."""
    sig = np.asarray(sigmag_amode, dtype=np.float64)
    alogsig = np.log(sig)
    return dict(
        alogsig=jnp.asarray(alogsig),
        exp45logsig=jnp.asarray(np.exp(4.5 * alogsig * alogsig)),
        f1=jnp.asarray(0.5 * np.exp(2.5 * alogsig * alogsig)),
        f2=jnp.asarray(1.0 + 0.25 * alogsig),
        aten=2.0 * MWH2O * SURFTEN / (RGAS * T0_ACT * RHOH2O),
    )


def _maxsat(zeta, eta, smc, f1, f2):
    """maxsat: maximum supersaturation over competing modes.
    zeta/eta/smc are (n, nmodes). Reproduces the Fortran scan: if ALL
    modes are weakly forced, smax = 1e-20; otherwise sum over all
    modes (eta <= 1e-20 lanes contribute 1e20)."""
    weak = (zeta > 1.0e5 * eta) | (smc * smc > 1.0e5 * eta)
    all_weak = jnp.all(weak, axis=1)
    eta_ok = eta > 1.0e-20
    g1 = zeta / W(eta_ok, eta, 1.0)
    g1 = jnp.sqrt(g1) * g1
    g2 = smc / jnp.sqrt(eta + 3.0 * zeta)
    g2 = jnp.sqrt(g2) * g2
    term = W(eta_ok, (f1[None, :] * g1 + f2[None, :] * g2)
             / (smc * smc), 1.0e20)
    # sequential mode accumulation (Fortran loop order, no pairwise)
    ssum = term[:, 0]
    for m in range(1, term.shape[1]):
        ssum = ssum + term[:, m]
    smax = 1.0 / jnp.sqrt(ssum)
    return W(all_weak, 1.0e-20, smax)


def actdrop_mam_calc(wbar, tair, rhoair, na, volume, hygro, in_cloud,
                     smax_f, ap):
    """actdrop_mam_calc with sigw = wdiab = 0 (zm_mphy configuration):
    the single-updraft branch. Returns fn (n, nmodes), the number
    fraction of aerosols activated (the only output zm_mphy uses).
    All inputs vectorized over the leading axis; masked-out lanes
    must carry finite dummies."""
    n, nmode = na.shape
    alogsig = ap["alogsig"]
    aten = ap["aten"]

    pres = RAIR * rhoair * tair
    diff0 = 0.211e-4 * (P0_ACT / pres) * _divs(tair, T0_ACT) ** 1.94
    conduct0 = (5.69 + 0.017 * (tair - T0_ACT)) * 4.186e2 * 1.0e-5
    r = _qsat_pa(tair, pres)
    qs = r["qs"]
    dqsdt = LATVAP / (RH2O * tair * tair) * qs
    alpha = GRAVIT * (LATVAP / (CPAIR * RH2O * tair * tair)
                      - 1.0 / (RAIR * tair))
    gamma_ = (1.0 + LATVAP / CPAIR * dqsdt) / (rhoair * qs)
    wmaxf = 10.0
    etafactor2max = 1.0e10 / (alpha * wmaxf) ** 1.5
    grow = 1.0 / (RHOH2O / (diff0 * rhoair * qs)
                  + LATVAP * RHOH2O / (conduct0 * tair)
                  * (LATVAP / (RH2O * tair) - 1.0))
    sqrtg = jnp.sqrt(grow)[:, None] * jnp.ones((1, nmode))

    ok = (volume > 1.0e-39) & (na > 1.0e-39)
    na_s = W(ok, na, 1.0)
    amcube = 3.0 * volume / (4.0 * PI * ap["exp45logsig"][None, :] * na_s)
    beta = 2.0 * PI * RHOH2O * grow * gamma_
    etafactor2 = W(ok, 1.0 / (na_s * beta[:, None] * sqrtg),
                   etafactor2max[:, None])
    hyg_ok = hygro > 1.0e-10
    smc = W(ok,
            W(hyg_ok, 2.0 * aten
              * jnp.sqrt(aten / (27.0 * W(hyg_ok, hygro, 1.0)
                                 * W(ok, amcube, 1.0))), 100.0),
            1.0)
    lnsm = jnp.log(smc)

    # single updraft: wnuc = wbar (+wdiab=0); no activation if wnuc<=0
    wnuc = wbar
    active = wnuc > 0.0
    w = wbar

    # in-cloud smax
    smax_f_pos = smax_f > 0.0
    smax_ic = W(smax_f_pos,
                alpha * w / (2.0 * PI * RHOH2O * grow * gamma_
                             * W(smax_f_pos, smax_f, 1.0)),
                1.0e-20)

    # cloud-base smax via maxsat
    alw = alpha * W(active, wnuc, 1.0)
    sqrtalw = jnp.sqrt(alw)
    etafactor1 = alw * sqrtalw
    eta = etafactor1[:, None] * etafactor2
    zeta = TWOTHIRD * sqrtalw[:, None] * aten / sqrtg
    smax_cb = _maxsat(zeta, eta, smc, ap["f1"], ap["f2"])

    smax = W(in_cloud, smax_ic, smax_cb)
    lnsmax = jnp.log(smax)
    x = _bdiv(TWOTHIRD * (lnsm - lnsmax[:, None]),
              (SQ2 * alogsig)[None, :])
    fn = 0.5 * (1.0 - _erf(x))
    return W(active[:, None], fn, 0.0)


# ---------------------------------------------------------------------------
# nucleate_ice_conv (Liu & Penner 2005)
# ---------------------------------------------------------------------------
def _hetero(tc, ww, ns):
    """hetero: immersion freezing of soot+dust [#/cm3]."""
    a11, a12 = 0.0263, -0.0185
    a21, a22 = 2.758, 1.3221
    b11, b12 = -0.008, -0.0468
    b21, b22 = -0.2667, -1.4588
    lns = jnp.log(ns)
    lnw = jnp.log(ww)
    b = (a11 + b11 * lns) * lnw + (a12 + b12 * lns)
    c = a21 + b21 * lns
    nis = jnp.exp(a22) * ns ** b22 * jnp.exp(b * tc) * ww ** c
    nis = jnp.minimum(nis, ns)
    return nis, jnp.zeros_like(nis)


def _hf(t, ww, rh, subgrid, na):
    """hf: homogeneous freezing of sulfate [#/cm3]."""
    a1_fast, a21_fast, a22_fast = 0.0231, -1.6387, -6.045
    b1_fast, b21_fast, b22_fast = -0.008, -0.042, -0.112
    c1_fast, c2_fast = 0.0739, 1.2372
    a1_slow, a2_slow = -0.3949, 1.282
    b1_slow, b2_slow, b3_slow = -0.0156, 0.0111, 0.0217
    c1_slow, c2_slow = 0.120, 2.312

    lnw = jnp.log(ww)
    a = 6.0e-4 * lnw + 6.6e-3
    b = 6.0e-2 * lnw + 1.052
    c = 1.68 * lnw + 129.35
    rhw = (a * t * t + b * t + c) * 0.01

    go = (t <= -37.0) & ((rh * subgrid) >= rhw)
    regm = 6.07 * lnw - 55.0
    fast = t >= regm
    warm64 = t > -64.0
    a2f = W(warm64, a21_fast, a22_fast)
    b2f = W(warm64, b21_fast, b22_fast)
    k1_fast = jnp.exp(a2f + b2f * t + c2_fast * lnw)
    k2_fast = a1_fast + b1_fast * t + c1_fast * lnw
    ni_fast = jnp.minimum(k1_fast * na ** k2_fast, na)
    k1_slow = jnp.exp(a2_slow + (b2_slow + b3_slow * lnw) * t
                      + c2_slow * lnw)
    k2_slow = a1_slow + b1_slow * t + c1_slow * lnw
    ni_slow = jnp.minimum(k1_slow * na ** k2_slow, na)
    return W(go, W(fast, ni_fast, ni_slow), 0.0)


def nucleati_conv(wbar, tair, relhum, cldn, qc, rhoair, so4_num,
                  dst_num, soot_num):
    """nucleati_conv with zm_microp = .true. (subgrid = 1). Returns
    nuci [#/kg] (the only output zm_mphy feeds back). Vectorized;
    masked lanes need finite inputs."""
    tc = tair - 273.15
    subgrid = 1.0

    outer = ((so4_num >= 1.0e-10) & ((soot_num + dst_num) >= 1.0e-10)
             & (cldn > 0.0))
    # safe values for log()s on masked lanes
    sn = W(outer, soot_num + dst_num, 1.0)
    wb = jnp.maximum(wbar, 1.0e-30)
    so4 = W(outer, jnp.maximum(so4_num, 1.0e-30), 1.0)

    cold = (tc <= -35.0) & ((relhum * svp_water(tair) / svp_ice(tair)
                             * subgrid) >= 1.2)
    a = -1.4938 * jnp.log(sn) + 12.884
    b = -10.41 * jnp.log(sn) - 67.69
    regm = a * jnp.log(wb) + b

    excl = (tc < -40.0) & (wbar > 1.0)
    nihf_only = _hf(tc, wb, relhum, subgrid, so4)
    nis, nid = _hetero(tc, wb, sn)

    # transition zone interpolation
    nihf_tr = _hf(regm - 5.0, wb, relhum, subgrid, so4)
    nis_tr, nid_tr = _hetero(regm, wb, sn)
    het_tr = nis_tr + nid_tr
    hf_le = nihf_tr <= het_tr
    ratio_safe = het_tr / W(nihf_tr > 0.0, nihf_tr, 1.0)
    n1_tr = W(hf_le, nihf_tr,
              het_tr * W(nihf_tr > 0.0, ratio_safe, 1.0)
              ** _divs(tc - regm, 5.0))

    hom_only = (tc > regm) & excl | (tc < regm - 5.0) \
        | ((tc >= regm - 5.0) & (tc <= regm) & excl)
    het_only = (tc > regm) & ~excl
    trans = ~hom_only & ~het_only

    n1 = W(hom_only, nihf_only, W(het_only, nis + nid, n1_tr))
    ni = W(outer & cold, n1, 0.0)

    # Meyers 1992 deposition/condensation nucleation in mixed clouds
    mixed = (tc < 0.0) & (tc > -37.0) & (qc > 1.0e-12)
    esl = svp_water(tair)
    esi = svp_ice(tair)
    deles = esl - esi
    nimey = W(mixed, 1.0e-3 * jnp.exp(12.96 * deles / esi - 0.639), 0.0)

    nuci = ni + nimey
    nuci = W((nuci > 9999.0) | (nuci < 0.0), 0.0, nuci)
    return nuci * 1.0e6 / rhoair


# ---------------------------------------------------------------------------
# zm_mphy
# ---------------------------------------------------------------------------
def zm_mphy(msg, jb, jt, jlcl, su, qu, mu, du, eu, zf, pm, te, qe,
            gamhat, eps0, cmel, cmei, aero, deltat, auto_fac,
            accr_fac, dcs, grav, cp, rd, lamc0=LAMBDADPCU,
            pgam0=MUCON, mp=None, ap=None):
    """zm_microphysics.F90 zm_mphy. All 2D inputs (ncol, pver)
    C-ordered, level 0 = top; zf (ncol, pver+1); jb/jt/jlcl 0-based;
    msg keeps its Fortran VALUE; pm in hPa. aero: dict with keys
    scheme='modal', nspec (list), mode_accum/mode_aitken/mode_coarse,
    coarse_dust/coarse_nacl/coarse_so4 (ALL 0-based), sigmag_aitken,
    specdens/spechygro (nspecmx, nmodes), voltonumblo/voltonumbhi
    (nmodes,), numg (ncol, pver, nmodes) [#/kg], mmrg (ncol, pver,
    nspecmx, nmodes) [kg/kg], dgnumg (ncol, pver, nmodes) [m],
    sigmag_amode (nmodes,). Returns dict of the 22 state fields and
    60 diagnostic fields (golden STATE_FIELDS/DIAG_FIELDS names)."""
    if aero["scheme"] != "modal":
        raise NotImplementedError(
            "zm_mphy port: modal aerosols only (EAMv3 production); the "
            "bulk ndrop_bam path is out of scope (PORTING_PLAN row 9)")
    if mp is None:
        mp = make_mphyi()
    if ap is None:
        ap = make_actdrop_params(aero["sigmag_amode"])

    su = jnp.asarray(su, dtype=jnp.float64)
    qu = jnp.asarray(qu, dtype=jnp.float64)
    mu = jnp.asarray(mu, dtype=jnp.float64)
    du = jnp.asarray(du, dtype=jnp.float64)
    eu = jnp.asarray(eu, dtype=jnp.float64)
    zf = jnp.asarray(zf, dtype=jnp.float64)
    pm = jnp.asarray(pm, dtype=jnp.float64)
    te = jnp.asarray(te, dtype=jnp.float64)
    qe = jnp.asarray(qe, dtype=jnp.float64)
    gamhat = jnp.asarray(gamhat, dtype=jnp.float64)
    eps0 = jnp.asarray(eps0, dtype=jnp.float64)
    cmel = jnp.asarray(cmel, dtype=jnp.float64)
    cmei = jnp.asarray(cmei, dtype=jnp.float64)
    jb = jnp.asarray(jb)
    jt = jnp.asarray(jt)
    jlcl = jnp.asarray(jlcl)

    n, pver = te.shape
    nmodes = int(aero["numg"].shape[2])
    numg = jnp.asarray(aero["numg"], dtype=jnp.float64)
    mmrg = jnp.asarray(aero["mmrg"], dtype=jnp.float64)
    dgnumg = jnp.asarray(aero["dgnumg"], dtype=jnp.float64)
    specdens = jnp.asarray(aero["specdens"], dtype=jnp.float64)
    spechygro = jnp.asarray(aero["spechygro"], dtype=jnp.float64)
    v2nlo = jnp.asarray(aero["voltonumblo"], dtype=jnp.float64)
    v2nhi = jnp.asarray(aero["voltonumbhi"], dtype=jnp.float64)
    nspec = list(aero["nspec"])
    m_acc = int(aero["mode_accum"])
    m_ait = int(aero["mode_aitken"])
    m_crs = int(aero["mode_coarse"])
    l_dst = int(aero["coarse_dust"])
    l_ncl = int(aero["coarse_nacl"])
    l_so4 = int(aero["coarse_so4"])
    sigmag_aitken = float(aero["sigmag_aitken"])

    xlf = mp["xlf"]
    ci, di = mp["ci"], mp["di"]
    cs, ds = mp["cs"], mp["ds"]
    cg, dg = mp["cg"], mp["dg"]
    rhosu = mp["rhosu"]
    mi0 = mp["mi0"]
    mmult = mp["mmult"]

    # parameters
    omsm = 0.99999
    zfacbuo = 0.5 / (1.0 + 0.5)
    cwdrag = 1.875 * 0.506
    retv = 0.608

    zeros = jnp.zeros((n, pver))

    # ---- initialization (Fortran 654-777) ----
    q = jnp.array(qu)
    tu = su - grav / cp * zf[:, :pver]
    t = su - grav / cp * zf[:, :pver]
    p = 100.0 * pm
    wu = jnp.zeros((n, pver))
    zkine = jnp.zeros((n, pver))
    arcf = jnp.zeros((n, pver))
    zbuo = jnp.zeros((n, pver))
    nc = jnp.zeros((n, pver))
    ni = jnp.zeros((n, pver))
    qc = jnp.zeros((n, pver))
    qi = jnp.zeros((n, pver))
    ncde = jnp.zeros((n, pver))
    nide = jnp.zeros((n, pver))
    nsde = jnp.zeros((n, pver))
    qcde = jnp.zeros((n, pver))
    qide = jnp.zeros((n, pver))
    qnide = jnp.zeros((n, pver))
    rprd = jnp.zeros((n, pver))
    sprd = jnp.zeros((n, pver))
    frz = jnp.zeros((n, pver))
    qr = jnp.zeros((n, pver))
    qni = jnp.zeros((n, pver))
    qg = jnp.zeros((n, pver))
    nr = jnp.zeros((n, pver))
    ns = jnp.zeros((n, pver))
    ng = jnp.zeros((n, pver))
    fhmrm = jnp.zeros((n, pver))
    fholm = jnp.zeros((n, pver))
    fholn = jnp.zeros((n, pver))
    # lamc/pgam are intent(inout): scalar seeds or full arrays
    # (persistent across the two zm_mphy calls of zm_cloud_properties)
    lamc = jnp.zeros((n, pver)) + jnp.asarray(lamc0, dtype=jnp.float64)
    pgam = jnp.zeros((n, pver)) + jnp.asarray(pgam0, dtype=jnp.float64)

    diag = {nm: jnp.zeros((n, pver)) for nm in _DIAG_NAMES}

    # time-varying parameters (Fortran 780-818)
    kidx = jnp.arange(pver)[None, :]
    p_km1 = jnp.concatenate([p[:, :1], p[:, :-1]], axis=1)
    te_km1 = jnp.concatenate([te[:, :1], te[:, :-1]], axis=1)
    qe_km1 = jnp.concatenate([qe[:, :1], qe[:, :-1]], axis=1)
    top = kidx == 0
    rhoh = W(top, p / (t * rd), 0.5 * (p + p_km1) / (t * rd))
    t_kp1 = jnp.concatenate([t[:, 1:], t[:, -1:]], axis=1)
    rhom = W(top, p / (t * rd),
             W(kidx == pver - 1, p / (rd * t),
               2.0 * p / (rd * (t + t_kp1))))
    th = W(top, te, 0.5 * (te + te_km1))
    qh = W(top, qe, 0.5 * (qe + qe_km1))
    dz = W(top, zf[:, :pver] - zf[:, 1:],
           jnp.concatenate([zf[:, :1], zf[:, :pver - 1]], axis=1)
           - zf[:, :pver])
    ph = W(top, p, 0.5 * (p + p_km1))
    dv = 8.794e-5 * t ** 1.81 / ph
    mua = 1.496e-6 * t ** 1.5 / (t + 120.0)
    rho = rhoh
    arn = AR * (rhosu / rho) ** 0.54
    asn = AS_ * (rhosu / rho) ** 0.54
    acn = AC * (rhosu / rho) ** 0.54
    ain = AI * (rhosu / rho) ** 0.54
    agn = AG * (rhosu / rho) ** 0.54

    # modal aerosol init (Fortran 820-837): ntaer accumulates over
    # modes sequentially (ntaer + numg*rhom per mode, Fortran order)
    ntaer = jnp.zeros((n, pver))
    for m in range(nmodes):
        ntaer = ntaer + numg[:, :, m] * rhom

    # ltrue (Fortran 876-881): qc/qi are all zero at this point
    ltrue = jnp.any((cmel >= QSMALL) | (cmei >= QSMALL), axis=1)

    kqc = jnp.zeros(n, dtype=jnp.int64)   # 0-based; Fortran init 1
    kqi = jnp.zeros(n, dtype=jnp.int64)
    lcbase = jnp.ones(n, dtype=bool)
    libase = jnp.ones(n, dtype=bool)

    # per-level work arrays persisted across levels and iterations
    lami = jnp.zeros((n, pver))
    n0i = jnp.zeros((n, pver))
    lams = jnp.zeros((n, pver))
    n0s = jnp.zeros((n, pver))
    lamg = jnp.zeros((n, pver))
    n0g = jnp.zeros((n, pver))
    lamr = jnp.zeros((n, pver))
    n0r = jnp.zeros((n, pver))
    cdist1 = jnp.zeros((n, pver))
    umr = jnp.zeros((n, pver))
    unr = jnp.zeros((n, pver))
    ums = jnp.zeros((n, pver))
    uns = jnp.zeros((n, pver))
    umg = jnp.zeros((n, pver))
    ung = jnp.zeros((n, pver))
    prf = jnp.zeros((n, pver))
    pnrf = jnp.zeros((n, pver))
    psf = jnp.zeros((n, pver))
    pnsf = jnp.zeros((n, pver))
    pgf = jnp.zeros((n, pver))
    pngf = jnp.zeros((n, pver))
    qctend = jnp.zeros((n, pver))
    qitend = jnp.zeros((n, pver))
    qnitend = jnp.zeros((n, pver))
    qrtend = jnp.zeros((n, pver))
    qgtend = jnp.zeros((n, pver))
    nctend = jnp.zeros((n, pver))
    nitend = jnp.zeros((n, pver))
    nrtend = jnp.zeros((n, pver))
    nstend = jnp.zeros((n, pver))
    ngtend = jnp.zeros((n, pver))
    qcic = jnp.zeros((n, pver))
    qiic = jnp.zeros((n, pver))
    qniic = jnp.zeros((n, pver))
    qric = jnp.zeros((n, pver))
    qgic = jnp.zeros((n, pver))
    ncic = jnp.zeros((n, pver))
    niic = jnp.zeros((n, pver))
    nsic = jnp.zeros((n, pver))
    nric = jnp.zeros((n, pver))
    ngic = jnp.zeros((n, pver))
    dum2l = jnp.zeros((n, pver))
    dum2i = jnp.zeros((n, pver))

    rows = jnp.arange(n)
    mtime = deltat / 900.0
    mtimec = deltat / 900.0

    # static species mask: contributions from l >= nspec(m) are zero
    spec_mask_np = np.zeros((specdens.shape[0], nmodes), dtype=bool)
    for m in range(nmodes):
        spec_mask_np[:nspec[m], m] = True
    spec_mask = jnp.asarray(spec_mask_np)

    for it in (1, 2):
        first = it == 1
        # sub-step re-init (Fortran 941-1051), ltrue columns only;
        # non-ltrue columns are all-zero anyway so a full reset of the
        # per-iteration arrays is exact
        qctend = zeros; qitend = zeros; qnitend = zeros  # noqa: E702
        qrtend = zeros; qgtend = zeros; nctend = zeros   # noqa: E702
        nitend = zeros; nrtend = zeros; nstend = zeros   # noqa: E702
        ngtend = zeros
        rprd = zeros; sprd = zeros; frz = zeros          # noqa: E702
        qniic = zeros; qric = zeros; qgic = zeros        # noqa: E702
        nsic = zeros; nric = zeros; ngic = zeros         # noqa: E702
        qiic = zeros; qcic = zeros; niic = zeros         # noqa: E702
        ncic = zeros
        dum2l = zeros; dum2i = zeros                     # noqa: E702
        fholm = zeros; fholn = zeros; fhmrm = zeros      # noqa: E702
        for nm in _DIAG_NAMES:
            diag[nm] = zeros
        ncadj = zeros
        niadj = zeros

        for kf in range(pver, msg + 1, -1):     # Fortran k
            c = kf - 1                          # python index of k
            cm = c - 1                          # python index of k-1
            act = (ltrue & (c > jt) & (c <= jb) & (eps0 > 0.0)
                   & (mu[:, c] > 0.0) & (mu[:, cm] > 0.0))
            if not bool(jnp.any(act)):
                # cleanups below the big if still run for every k
                (qni, ns, qr, nr, qg, ng, qi, ni, qc, nc) = _cleanup(
                    ltrue, cm, qni, ns, qr, nr, qg, ng, qi, ni, qc, nc)
                continue

            dzc = dz[:, c]
            muc = mu[:, c]
            mucm = mu[:, cm]

            if first:
                umr = umr.at[:, c].set(W(act, 0.0, umr[:, c]))
                unr = unr.at[:, c].set(W(act, 0.0, unr[:, c]))
                ums = ums.at[:, c].set(W(act, 0.0, ums[:, c]))
                uns = uns.at[:, c].set(W(act, 0.0, uns[:, c]))
                umg = umg.at[:, c].set(W(act, 0.0, umg[:, c]))
                ung = ung.at[:, c].set(W(act, 0.0, ung[:, c]))
                prf = prf.at[:, c].set(W(act, 0.0, prf[:, c]))
                pnrf = pnrf.at[:, c].set(W(act, 0.0, pnrf[:, c]))
                psf = psf.at[:, c].set(W(act, 0.0, psf[:, c]))
                pnsf = pnsf.at[:, c].set(W(act, 0.0, pnsf[:, c]))
                pgf = pgf.at[:, c].set(W(act, 0.0, pgf[:, c]))
                pngf = pngf.at[:, c].set(W(act, 0.0, pngf[:, c]))

            # ---- in-updraft values (Fortran 1084-1197) ----
            if first:
                qcic = qcic.at[:, c].set(W(act, qc[:, c], qcic[:, c]))
                qiic = qiic.at[:, c].set(W(act, qi[:, c], qiic[:, c]))
                ncic = ncic.at[:, c].set(W(act, nc[:, c], ncic[:, c]))
                niic = niic.at[:, c].set(W(act, ni[:, c], niic[:, c]))
                qniic = qniic.at[:, c].set(W(act, qni[:, c], qniic[:, c]))
                qric = qric.at[:, c].set(W(act, qr[:, c], qric[:, c]))
                nsic = nsic.at[:, c].set(W(act, ns[:, c], nsic[:, c]))
                nric = nric.at[:, c].set(W(act, nr[:, c], nric[:, c]))
                qgic = qgic.at[:, c].set(W(act, qg[:, c], qgic[:, c]))
                ngic = ngic.at[:, c].set(W(act, ng[:, c], ngic[:, c]))
            else:
                mc = act & (c <= kqc)
                qcic = qcic.at[:, c].set(W(mc, qc[:, c], qcic[:, c]))
                ncic = ncic.at[:, c].set(W(mc, nc[:, c], ncic[:, c]))
                # rain falling from above: Fortran kk = k .. jt+3
                # descending (1-based); lf is Fortran kk, so with
                # 0-based jt the bound is lf >= jt+4 and array reads
                # sit at python index lf-2 (= kk-1)
                flxrm = jnp.zeros(n); mvtrm = jnp.zeros(n)  # noqa: E702
                flxrn = jnp.zeros(n); mvtrn = jnp.zeros(n)  # noqa: E702
                jtmin = int(jnp.min(jt))
                for lf in range(kf, jtmin + 3, -1):
                    lm = lf - 2                 # python index kk-1
                    mk = mc & (lf >= jt + 4) & (qr[:, lm] > 0.0)
                    flxrm = flxrm + W(mk, umr[:, lm] * qr[:, lm]
                                      * arcf[:, lm], 0.0)
                    flxrn = flxrn + W(mk, unr[:, lm] * nr[:, lm]
                                      * arcf[:, lm], 0.0)
                    mvtrm = mvtrm + W(mk, umr[:, lm] * arcf[:, lm], 0.0)
                    mvtrn = mvtrn + W(mk, unr[:, lm] * arcf[:, lm], 0.0)
                qric = qric.at[:, c].set(
                    W(mc, W(mvtrm > 0,
                            (qr[:, c] * muc + flxrm)
                            / W(mvtrm > 0, muc + mvtrm, 1.0),
                            qr[:, c]), qric[:, c]))
                nric = nric.at[:, c].set(
                    W(mc, W(mvtrn > 0,
                            (nr[:, c] * muc + flxrn)
                            / W(mvtrn > 0, muc + mvtrn, 1.0),
                            nr[:, c]), nric[:, c]))
                meq = act & (c == kqc)
                qcic = qcic.at[:, c].set(W(meq, qc[:, cm], qcic[:, c]))
                ncic = ncic.at[:, c].set(W(meq, nc[:, cm], ncic[:, c]))

                mi_ = act & (c <= kqi)
                qiic = qiic.at[:, c].set(W(mi_, qi[:, c], qiic[:, c]))
                niic = niic.at[:, c].set(W(mi_, ni[:, c], niic[:, c]))
                # snow falling from above
                flxsm = jnp.zeros(n); mvtsm = jnp.zeros(n)  # noqa: E702
                flxsn = jnp.zeros(n); mvtsn = jnp.zeros(n)  # noqa: E702
                for lf in range(kf, jtmin + 3, -1):
                    lm = lf - 2
                    mk = mi_ & (lf >= jt + 4) & (qni[:, lm] > 0.0)
                    flxsm = flxsm + W(mk, ums[:, lm] * qni[:, lm]
                                      * arcf[:, lm], 0.0)
                    mvtsm = mvtsm + W(mk, ums[:, lm] * arcf[:, lm], 0.0)
                    flxsn = flxsn + W(mk, uns[:, lm] * ns[:, lm]
                                      * arcf[:, lm], 0.0)
                    mvtsn = mvtsn + W(mk, uns[:, lm] * arcf[:, lm], 0.0)
                qniic = qniic.at[:, c].set(
                    W(mi_, W(mvtsm > 0,
                             (qni[:, c] * muc + flxsm)
                             / W(mvtsm > 0, muc + mvtsm, 1.0),
                             qni[:, c]), qniic[:, c]))
                nsic = nsic.at[:, c].set(
                    W(mi_, W(mvtsn > 0,
                             (ns[:, c] * muc + flxsn)
                             / W(mvtsn > 0, muc + mvtsn, 1.0),
                             ns[:, c]), nsic[:, c]))
                # graupel falling from above
                flxgm = jnp.zeros(n); mvtgm = jnp.zeros(n)  # noqa: E702
                flxgn = jnp.zeros(n); mvtgn = jnp.zeros(n)  # noqa: E702
                for lf in range(kf, jtmin + 3, -1):
                    lm = lf - 2
                    mk = mi_ & (lf >= jt + 4) & (qg[:, lm] > 0.0)
                    flxgm = flxgm + W(mk, umg[:, lm] * qg[:, lm]
                                      * arcf[:, lm], 0.0)
                    mvtgm = mvtgm + W(mk, umg[:, lm] * arcf[:, lm], 0.0)
                    flxgn = flxgn + W(mk, ung[:, lm] * ng[:, lm]
                                      * arcf[:, lm], 0.0)
                    mvtgn = mvtgn + W(mk, ung[:, lm] * arcf[:, lm], 0.0)
                qgic = qgic.at[:, c].set(
                    W(mi_, W(mvtgm > 0,
                             (qg[:, c] * muc + flxgm)
                             / W(mvtgm > 0, muc + mvtgm, 1.0),
                             qg[:, c]), qgic[:, c]))
                ngic = ngic.at[:, c].set(
                    W(mi_, W(mvtgn > 0,
                             (ng[:, c] * muc + flxgn)
                             / W(mvtgn > 0, muc + mvtgn, 1.0),
                             ng[:, c]), ngic[:, c]))
                meqi = act & (c == kqi)
                qiic = qiic.at[:, c].set(W(meqi, qi[:, cm], qiic[:, c]))
                niic = niic.at[:, c].set(W(meqi, ni[:, cm], niic[:, c]))

            # ---- boundary conditions (Fortran 1204-1225, it=1) ----
            if first:
                bc_c = act & (cmel[:, cm] > QSMALL) & lcbase
                kqc = W(bc_c, c, kqc)
                lcbase = lcbase & ~bc_c
                den_c = mucm + dzc * du[:, cm]
                qcic_bc = dzc * cmel[:, cm] / W(bc_c, den_c, 1.0)
                qcic = qcic.at[:, c].set(W(bc_c, qcic_bc, qcic[:, c]))
                ncic = ncic.at[:, c].set(
                    W(bc_c, _divs(qcic[:, c],
                                  4.0 / 3.0 * PI * 25.0e-6 ** 3 * RHOW),
                      ncic[:, c]))

                bc_i1 = act & (qiic[:, c] > QSMALL) & libase
                bc_i2 = (act & ~bc_i1 & (cmei[:, cm] > QSMALL)
                         & (cmei[:, c] < QSMALL) & (c <= jb) & libase)
                kqi = W(bc_i1 | bc_i2, c, kqi)
                libase = libase & ~(bc_i1 | bc_i2)
                qiic_bc = dzc * cmei[:, cm] / W(bc_i2, den_c, 1.0)
                qiic = qiic.at[:, c].set(W(bc_i2, qiic_bc, qiic[:, c]))
                niic = niic.at[:, c].set(
                    W(bc_i2, _divs(qiic[:, c],
                                   4.0 / 3.0 * PI * 15.0e-6 ** 3 * RHOI),
                      niic[:, c]))

            # ---- cloud ice/water size distributions (1232-1298) ----
            mi_sd = act & (qiic[:, c] >= QSMALL)
            niic = niic.at[:, c].set(
                W(mi_sd, jnp.minimum(niic[:, c], qiic[:, c] * 1.0e20),
                  niic[:, c]))
            qiic_s = W(mi_sd, qiic[:, c], 1.0)
            lami_c = (mp["g_1di"] * ci
                      * W(mi_sd, niic[:, c], 1.0) / qiic_s) ** (1.0 / di)
            lammax_i = 1.0 / 10.0e-6
            lammin_i = 1.0 / (2.0 * dcs)
            lo = mi_sd & (lami_c < lammin_i)
            hi = mi_sd & (lami_c > lammax_i)
            lami_c = W(lo, lammin_i, W(hi, lammax_i, lami_c))
            n0i_c = W(lo | hi,
                      _divs(lami_c ** (di + 1.0) * qiic_s,
                            ci * mp["g_1di"]),
                      W(mi_sd, niic[:, c], 0.0) * lami_c)
            niic = niic.at[:, c].set(
                W(lo | hi, n0i_c / lami_c, niic[:, c]))
            lami = lami.at[:, c].set(W(mi_sd, lami_c, W(act, 0.0,
                                                        lami[:, c])))
            n0i = n0i.at[:, c].set(W(mi_sd, n0i_c, W(act, 0.0,
                                                     n0i[:, c])))

            mc_sd = act & (qcic[:, c] >= QSMALL)
            ncic = ncic.at[:, c].set(
                W(mc_sd, jnp.minimum(ncic[:, c], qcic[:, c] * 1.0e20),
                  ncic[:, c]))
            pg = 0.0005714 * (_divs(W(mc_sd, ncic[:, c], 0.0), 1.0e6)
                              * rho[:, c]) + 0.2714
            pg = 1.0 / (pg ** 2) - 1.0
            pg = jnp.clip(pg, 2.0, 15.0)
            pgam = pgam.at[:, c].set(W(mc_sd, pg, pgam[:, c]))
            qcic_s = W(mc_sd, qcic[:, c], 1.0)
            g_pg4 = _gamma(pg + 4.0)
            g_pg1 = _gamma(pg + 1.0)
            lamc_c = (PI / 6.0 * RHOW * W(mc_sd, ncic[:, c], 1.0)
                      * g_pg4 / (qcic_s * g_pg1)) ** (1.0 / 3.0)
            lammin_c = _divs(pg + 1.0, 40.0e-6)
            lammax_c = _divs(pg + 1.0, 1.0e-6)
            lo = mc_sd & (lamc_c < lammin_c)
            hi = mc_sd & (lamc_c > lammax_c)
            lamc_c = W(lo, lammin_c, W(hi, lammax_c, lamc_c))
            ncic = ncic.at[:, c].set(
                W(lo | hi, 6.0 * lamc_c ** 3 * qcic_s * g_pg1
                  / (PI * RHOW * g_pg4), ncic[:, c]))
            lamc = lamc.at[:, c].set(W(mc_sd, lamc_c,
                                       W(act, 0.0, lamc[:, c])))
            cdist1 = cdist1.at[:, c].set(
                W(mc_sd, ncic[:, c] / g_pg1, W(act, 0.0, cdist1[:, c])))

            # boundary zeroing (1300-1310)
            mb = act & (kqc == c)
            qc = qc.at[:, c].set(W(mb, 0.0, qc[:, c]))
            nc = nc.at[:, c].set(W(mb, 0.0, nc[:, c]))
            mb = act & (kqi == c)
            qi = qi.at[:, c].set(W(mb, 0.0, qi[:, c]))
            ni = ni.at[:, c].set(W(mb, 0.0, ni[:, c]))

            # ---- autoconversion of liquid (1321-1348) ----
            mauto = act & (qcic[:, c] >= 1.0e-8)
            qcic_s = W(mauto, qcic[:, c], 1.0)
            ncic_s = W(mauto, ncic[:, c], 1.0)
            prc = W(mauto, auto_fac * 30500.0 * qcic_s ** 3.19
                    * (_divs(ncic_s, 1.0e6) * rho[:, c]) ** (-1.2), 0.0)
            nprc1 = W(mauto, prc / (qcic_s / ncic_s), 0.0)
            nprc = W(mauto, prc * (1.0 / mp["droplet_mass_25um"]), 0.0)

            if first:
                mb = act & (c == kqc)
                qric = qric.at[:, c].set(W(mb, _divs(prc * dzc, 0.55),
                                           qric[:, c]))
                nric = nric.at[:, c].set(W(mb, _divs(nprc * dzc, 0.55),
                                           nric[:, c]))
                qr = qr.at[:, c].set(W(mb, 0.0, qr[:, c]))
                nr = nr.at[:, c].set(W(mb, 0.0, nr[:, c]))

            # ---- autoconversion of ice to snow (1354-1374) ----
            mia = act & (t[:, c] <= 273.15) & (qiic[:, c] >= QSMALL)
            lami_s = W(mia, lami[:, c], 1.0)
            nprci = W(mia, n0i[:, c] / (lami_s * 180.0)
                      * jnp.exp(-lami_s * dcs), 0.0)
            prci = W(mia, _divs(PI * RHOI * n0i[:, c], 6.0 * 180.0)
                     * (dcs ** 3 / lami_s + 3.0 * dcs ** 2 / lami_s ** 2
                        + 6.0 * dcs / lami_s ** 3 + 6.0 / lami_s ** 4)
                     * jnp.exp(-lami_s * dcs), 0.0)
            if first:
                mb = act & (c == kqi)
                qniic = qniic.at[:, c].set(W(mb, prci * dzc * 0.25,
                                             qniic[:, c]))
                nsic = nsic.at[:, c].set(W(mb, nprci * dzc * 0.25,
                                           nsic[:, c]))
                qni = qni.at[:, c].set(W(mb, 0.0, qni[:, c]))
                ns = ns.at[:, c].set(W(mb, 0.0, ns[:, c]))

            # zero number when mixing ratio is zero (1377-1394)
            mz = act & (qniic[:, c] < QSMALL)
            qniic = qniic.at[:, c].set(W(mz, 0.0, qniic[:, c]))
            nsic = nsic.at[:, c].set(W(mz, 0.0, nsic[:, c]))
            mz = act & (qric[:, c] < QSMALL)
            qric = qric.at[:, c].set(W(mz, 0.0, qric[:, c]))
            nric = nric.at[:, c].set(W(mz, 0.0, nric[:, c]))
            mz = act & (qgic[:, c] < QSMALL)
            qgic = qgic.at[:, c].set(W(mz, 0.0, qgic[:, c]))
            ngic = ngic.at[:, c].set(W(mz, 0.0, ngic[:, c]))
            nric = nric.at[:, c].set(
                W(act, jnp.maximum(nric[:, c], 0.0), nric[:, c]))
            nsic = nsic.at[:, c].set(
                W(act, jnp.maximum(nsic[:, c], 0.0), nsic[:, c]))
            ngic = ngic.at[:, c].set(
                W(act, jnp.maximum(ngic[:, c], 0.0), ngic[:, c]))

            # ---- precip size distributions (1400-1492) ----
            (lamr, n0r, nric, umr, unr) = _rain_sd(
                act, c, qric, nric, lamr, n0r, umr, unr, arn, rho, mp,
                pre=True)
            (lams, n0s, nsic, ums, uns) = _snow_sd(
                act, c, c, qniic, nsic, lams, n0s, ums, uns, asn, rho,
                mp, dcs, sp054=False)
            (lamg, n0g, ngic, umg, ung) = _graupel_sd(
                act, c, c, qgic, ngic, lamg, n0g, umg, ung, agn, rho,
                mp, sp054=False)

            # ---- snow self-aggregation (1500-1508) ----
            msa = act & (qniic[:, c] >= QSMALL) & (t[:, c] <= 273.15)
            nsagg = W(msa, -1108.0 * asn[:, c] * EII
                      * PI ** ((1.0 - BS) / 3.0)
                      * RHOSN ** ((-2.0 - BS) / 3.0)
                      * rho[:, c] ** ((2.0 + BS) / 3.0)
                      * W(msa, qniic[:, c], 1.0) ** ((2.0 + BS) / 3.0)
                      * (W(msa, nsic[:, c], 1.0)
                         * rho[:, c]) ** ((4.0 - BS) / 3.0)
                      / (4.0 * 720.0 * rho[:, c]), 0.0)

            # ---- accretion of droplets by snow + HM (1518-1559) ----
            mps = (act & (qniic[:, c] >= QSMALL) & (t[:, c] <= 273.15)
                   & (qcic[:, c] >= QSMALL))
            lamc_s = W(mps, lamc[:, c], 1.0)
            lams_s = W(mps, lams[:, c], 1.0)
            dc0 = (pgam[:, c] + 1.0) / lamc_s
            ds0 = 1.0 / lams_s
            dum = dc0 * dc0 * uns[:, c] * RHOW / (9.0 * mua[:, c] * ds0)
            eci = jnp.clip(dum * dum / ((dum + 0.4) * (dum + 0.4)),
                           0.0, 1.0)
            psacws = W(mps, PI / 4.0 * asn[:, c] * qcic[:, c]
                       * rho[:, c] * n0s[:, c] * eci * mp["g_bs3"]
                       / lams_s ** (BS + 3.0), 0.0)
            npsacws = W(mps, PI / 4.0 * asn[:, c] * ncic[:, c]
                        * rho[:, c] * n0s[:, c] * eci * mp["g_bs3"]
                        / lams_s ** (BS + 3.0), 0.0)
            hm1 = act & (t[:, c] < 270.16) & (t[:, c] >= 268.16)
            hm2 = act & (t[:, c] < 268.16) & (t[:, c] >= 265.16)
            ni_secp = W(hm1, 3.5e8 * (270.16 - t[:, c]) / 2.0 * psacws,
                        W(hm2, _divs(3.5e8 * (t[:, c] - 265.16), 3.0)
                          * psacws, 0.0))
            nsacwi = ni_secp
            msacwi = W(hm1 | hm2, jnp.minimum(ni_secp * mi0, psacws),
                       0.0)
            psacws = jnp.maximum(0.0, psacws - ni_secp * mi0)

            # ---- riming conversion to graupel (1567-1588) ----
            mgs = (act & (psacws > 0.0) & (qniic[:, c] >= 0.1e-3)
                   & (qcic[:, c] >= 0.5e-3))
            dt_g = dzc / W(mgs, ums[:, c], 1.0)
            pgsacw = W(mgs, jnp.minimum(
                psacws, mp["cons17"] * dt_g * n0s[:, c] * qcic[:, c]
                * qcic[:, c] * asn[:, c] * asn[:, c] * rho[:, c]
                / (lams_s ** (2.0 * BS + 2.0))), 0.0)
            dum = jnp.maximum(RHOSN / (RHOG - RHOSN) * pgsacw, 0.0)
            nscng = W(mgs, _divs(dum, MG0), 0.0)
            psacws = W(mgs, psacws - pgsacw, psacws)

            # ---- cloud ice collecting droplets (1595-1612) ----
            mii = (act & (qiic[:, c] >= 1.0e-8) & (qcic[:, c] >= QSMALL)
                   & (1.0 / W(act & (lami[:, c] > 0), lami[:, c], 1.0)
                      >= 100.0e-6))
            lami_s2 = W(mii, lami[:, c], 1.0)
            psacwi = W(mii, mp["cons16"] * ain[:, c] * qcic[:, c]
                       * rho[:, c] * n0i[:, c] / lami_s2 ** (BI + 3.0),
                       0.0)
            npsacwi = W(mii, mp["cons16"] * ain[:, c] * ncic[:, c]
                        * rho[:, c] * n0i[:, c] / lami_s2 ** (BI + 3.0),
                        0.0)

            # ---- collection of droplets by graupel (1618-1629) ----
            mwg = (act & (qgic[:, c] >= 1.0e-8) & (qcic[:, c] >= QSMALL)
                   & (t[:, c] <= 273.15))
            lamg_s = W(mwg, lamg[:, c], 1.0)
            psacwg = W(mwg, mp["cons14"] * agn[:, c] * qcic[:, c]
                       * rho[:, c] * n0g[:, c] / lamg_s ** (BG + 3.0),
                       0.0)
            npsacwg = W(mwg, mp["cons14"] * agn[:, c] * ncic[:, c]
                        * rho[:, c] * n0g[:, c] / lamg_s ** (BG + 3.0),
                        0.0)

            # ---- accretion of rain by snow (1636-1696) ----
            mrs = (act & (qric[:, c] >= 1.0e-8)
                   & (qniic[:, c] >= 1.0e-8) & (t[:, c] <= 273.15))
            lamr_s = W(mrs, lamr[:, c], 1.0)
            lams_s3 = W(mrs, lams[:, c], 1.0)
            pracs = W(mrs, PI * PI * ECR
                      * (((1.2 * umr[:, c] - 0.95 * ums[:, c]) ** 2
                          + 0.08 * ums[:, c] * umr[:, c]) ** 0.5
                         * RHOW * rho[:, c] * n0r[:, c] * n0s[:, c]
                         * (5.0 / (lamr_s ** 6 * lams_s3)
                            + 2.0 / (lamr_s ** 5 * lams_s3 ** 2.0)
                            + 0.5 / (lamr_s ** 4 * lams_s3 ** 3))),
                      0.0)
            npracs = W(mrs, PI / 2.0 * rho[:, c] * ECR
                       * (1.7 * (unr[:, c] - uns[:, c]) ** 2.0
                          + 0.3 * unr[:, c] * uns[:, c]) ** 0.5
                       * n0r[:, c] * n0s[:, c]
                       * (1.0 / (lamr_s ** 3.0 * lams_s3)
                          + 1.0 / (lamr_s ** 2.0 * lams_s3 ** 2.0)
                          + 1.0 / (lamr_s * lams_s3 ** 3.0)), 0.0)
            mps2 = mrs & (qniic[:, c] >= 0.1e-3) & (qric[:, c] >= 0.1e-3)
            psacr = W(mps2, mp["cons31"]
                      * (((1.2 * umr[:, c] - 0.95 * ums[:, c]) ** 2
                          + 0.08 * ums[:, c] * umr[:, c]) ** 0.5
                         * rho[:, c] * n0r[:, c] * n0s[:, c]
                         / lams_s3 ** 3.0
                         * (5.0 / (lams_s3 ** 3.0 * lamr_s)
                            + 2.0 / (lams_s3 ** 2.0 * lamr_s ** 2.0)
                            + 0.5 / (lams_s3 * lamr_s ** 3.0))), 0.0)
            mgr = (act & (pracs > 0.0) & (qniic[:, c] >= 0.1e-3)
                   & (qric[:, c] >= 0.1e-3))
            dum = (mp["cons18"] * (4.0 / lams_s3) ** 3.0
                   * (4.0 / lams_s3) ** 3.0
                   / (mp["cons18"] * (4.0 / lams_s3) ** 3.0
                      * (4.0 / lams_s3) ** 3.0
                      + mp["cons19"] * (4.0 / lamr_s) ** 3.0
                      * (4.0 / lamr_s) ** 3.0))
            dum = jnp.clip(dum, 0.0, 1.0)
            pgracs = W(mgr, (1.0 - dum) * pracs, 0.0)
            ngracs = W(mgr, (1.0 - dum) * npracs, 0.0)
            pracs = W(mgr, pracs - pgracs, pracs)
            npracs = W(mgr, npracs - ngracs, npracs)
            psacr = W(mgr, psacr * (1.0 - dum), psacr)

            # ---- collection of rain by graupel (1702-1720) ----
            mrg = (act & (qric[:, c] >= _THR_1EM8_SP)
                   & (qgic[:, c] >= _THR_1EM8_SP))
            lamr_s2 = W(mrg, lamr[:, c], 1.0)
            lamg_s2 = W(mrg, lamg[:, c], 1.0)
            pracg = W(mrg, mp["cons41"]
                      * (((1.2 * umr[:, c] - 0.95 * umg[:, c]) ** 2.0
                          + 0.08 * umg[:, c] * umr[:, c]) ** 0.5
                         * rho[:, c] * n0r[:, c] * n0g[:, c]
                         / lamr_s2 ** 3
                         * (5.0 / (lamr_s2 ** 3.0 * lamg_s2)
                            + 2.0 / (lamr_s2 ** 2.0 * lamg_s2 ** 2.0)
                            + 0.5 / (lamr_s2 * lamg_s2 ** 3.0))), 0.0)
            npracg = W(mrg, mp["cons32"] * rho[:, c]
                       * (1.7 * (unr[:, c] - ung[:, c]) ** 2.0
                          + 0.3 * unr[:, c] * ung[:, c]) ** 0.5
                       * n0r[:, c] * n0g[:, c]
                       * (1.0 / (lamr_s2 ** 3.0 * lamg_s2)
                          + 1.0 / (lamr_s2 ** 2 * lamg_s2 ** 2.0)
                          + 1.0 / (lamr_s2 * lamg_s2 ** 3.0)), 0.0)

            # ---- rime-splintering, graupel HM (1732-1775) ----
            hm_win = (act & (t[:, c] < 270.16) & (t[:, c] >= 265.16)
                      & (qgic[:, c] >= 0.1e-3)
                      & ((qcic[:, c] >= 0.5e-3)
                         | (qric[:, c] >= 0.1e-3))
                      & ((psacwg > 0.0) | (pracg > 0.0)))
            fmult = W(t[:, c] >= 268.16, (270.16 - t[:, c]) / 2.0,
                      _divs(t[:, c] - 265.16, 3.0))
            m_w = hm_win & (psacwg > 0.0)
            nmultg = W(m_w, 35.0e4 * psacwg * fmult, 0.0)
            qmultg = W(m_w, jnp.minimum(nmultg * mmult, psacwg), 0.0)
            psacwg = W(m_w, psacwg - qmultg, psacwg)
            m_r = hm_win & (pracg > 0.0)
            nmultrg = W(m_r, 35.0e4 * pracg * fmult, 0.0)
            qmultrg = W(m_r, jnp.minimum(nmultrg * mmult, pracg), 0.0)
            pracg = W(m_r, pracg - qmultrg, pracg)

            # ---- heterogeneous rain freezing (1781-1792) ----
            mfr = act & (t[:, c] < 269.15) & (qric[:, c] >= QSMALL)
            lamr_s3 = W(mfr, lamr[:, c], 1.0)
            mnuccr = W(mfr, 20.0 * PI * PI * RHOW * nric[:, c] * BIMM
                       * (jnp.exp(AIMM * (273.15 - t[:, c])) - 1.0)
                       / lamr_s3 ** 3.0 / lamr_s3 ** 3.0, 0.0)
            nnuccr = W(mfr, PI * nric[:, c] * BIMM
                       * (jnp.exp(AIMM * (273.15 - t[:, c])) - 1.0)
                       / lamr_s3 ** 3.0, 0.0)

            # ---- accretion of liquid by rain (1799-1805) ----
            mar = act & (qric[:, c] >= QSMALL) & (qcic[:, c] >= QSMALL)
            pra = W(mar, accr_fac * 67.0
                    * (qcic[:, c] * qric[:, c]) ** 1.15, 0.0)
            # pra/(qcic/ncic): IEEE x/0=inf, y/inf=0 when ncic=0
            npra = W(mar, pra / (W(mar, qcic[:, c], 1.0)
                                 / W(mar, ncic[:, c], 1.0)), 0.0)

            # ---- rain self-collection (1811-1815) ----
            msc = act & (qric[:, c] >= QSMALL)
            nragg = W(msc, -8.0 * nric[:, c] * qric[:, c] * rho[:, c],
                      0.0)

            # ---- accretion of ice by snow (1822-1833) ----
            mis = (act & (qniic[:, c] >= QSMALL)
                   & (qiic[:, c] >= QSMALL) & (t[:, c] <= 273.15))
            lams_s4 = W(mis, lams[:, c], 1.0)
            prai = W(mis, PI / 4.0 * asn[:, c] * qiic[:, c] * rho[:, c]
                     * n0s[:, c] * EII * mp["g_bs3"]
                     / lams_s4 ** (BS + 3.0), 0.0)
            nprai = W(mis, PI / 4.0 * asn[:, c] * niic[:, c] * rho[:, c]
                      * n0s[:, c] * EII * mp["g_bs3"]
                      / lams_s4 ** (BS + 3.0), 0.0)

            # ---- rain-ice collisions (1840-1867) ----
            mri = (act & (qric[:, c] >= 1.0e-8)
                   & (qiic[:, c] >= 1.0e-8) & (t[:, c] <= 273.15))
            mri_g = mri & (qric[:, c] >= 0.1e-3)
            mri_s = mri & ~mri_g
            lamr_s4 = W(mri, lamr[:, c], 1.0)
            base_n = (mp["cons24"] * niic[:, c] * n0r[:, c] * arn[:, c]
                      / lamr_s4 ** (BR + 3.0) * rho[:, c])
            base_p = (mp["cons25"] * niic[:, c] * n0r[:, c] * arn[:, c]
                      / lamr_s4 ** (BR + 3.0) / lamr_s4 ** 3 * rho[:, c])
            base_q = (mp["cons24"] * qiic[:, c] * n0r[:, c] * arn[:, c]
                      / lamr_s4 ** (BR + 3.0) * rho[:, c])
            niacr = W(mri_g, base_n, 0.0)
            piacr = W(mri_g, base_p, 0.0)
            praci = W(mri_g, base_q, 0.0)
            niacrs = W(mri_s, base_n, 0.0)
            piacrs = W(mri_s, base_p, 0.0)
            pracis = W(mri_s, base_q, 0.0)

            # ---- fallout terms (1871-1876) ----
            prf = prf.at[:, c].set(W(act, -umr[:, c] * qric[:, c] / dzc,
                                     prf[:, c]))
            pnrf = pnrf.at[:, c].set(W(act, -unr[:, c] * nric[:, c]
                                       / dzc, pnrf[:, c]))
            psf = psf.at[:, c].set(W(act, -ums[:, c] * qniic[:, c]
                                     / dzc, psf[:, c]))
            pnsf = pnsf.at[:, c].set(W(act, -uns[:, c] * nsic[:, c]
                                       / dzc, pnsf[:, c]))
            pgf = pgf.at[:, c].set(W(act, -umg[:, c] * qgic[:, c]
                                     / dzc, pgf[:, c]))
            pngf = pngf.at[:, c].set(W(act, -ung[:, c] * ngic[:, c]
                                       / dzc, pngf[:, c]))

            # ---- vertical velocity (1881-1913) ----
            mjb = act & (c == jb)
            zkine = zkine.at[:, c].set(W(mjb, 0.5, zkine[:, c]))
            wu = wu.at[:, c].set(W(mjb, 1.0, wu[:, c]))
            zbuo_jb = ((tu[:, c] * (1.0 + retv * qu[:, c])
                        - th[:, c] * (1.0 + retv * qh[:, c]))
                       / (th[:, c] * (1.0 + retv * qh[:, c])))
            zbuo = zbuo.at[:, c].set(W(mjb, zbuo_jb, zbuo[:, c]))
            mnb = act & (c != jb)
            cp1 = min(c + 1, pver - 1)
            zbc = tu[:, c] * (1.0 + retv * qu[:, c] - qr[:, c]
                              - qni[:, c] - qi[:, c] - qc[:, c])
            zbe = th[:, c] * (1.0 + retv * qh[:, c])
            zbuo_c = (zbc - zbe) / zbe
            zbuo = zbuo.at[:, c].set(W(mnb, zbuo_c, zbuo[:, c]))
            zbuoc = (zbuo[:, c] + zbuo[:, cp1]) * 0.5
            zdkbuo = dz[:, cp1] * grav * zfacbuo * zbuoc
            zdken = jnp.minimum(
                0.99, (1.0 + cwdrag) * jnp.maximum(du[:, c], eu[:, c])
                * dz[:, cp1] / jnp.maximum(1.0e-10, mu[:, cp1]))
            zkine = zkine.at[:, c].set(
                W(mnb, (zkine[:, cp1] * (1.0 - zdken) + zdkbuo)
                  / (1.0 + zdken), zkine[:, c]))
            wu = wu.at[:, c].set(
                W(mnb, jnp.minimum(
                    15.0, jnp.sqrt(2.0 * jnp.maximum(0.1, zkine[:, c]))),
                  wu[:, c]))
            arcf = arcf.at[:, c].set(
                W(act, muc / W(act, wu[:, c], 1.0), arcf[:, c]))

            # ---- droplet activation (1920-2012) ----
            ntaerh_c = 0.5 * (ntaer[:, c] + ntaer[:, cm])
            mact = act & (qcic[:, c] >= QSMALL)
            vol = _bdiv(jnp.maximum(
                0.5 * (mmrg[:, c, :, :] + mmrg[:, cm, :, :]), 0.0),
                specdens[None, :, :])
            # zero contributions from species l >= nspec(m); accumulate
            # sequentially over species (Fortran l-loop order)
            vol = W(spec_mask[None], vol, 0.0)
            vaerosol = jnp.zeros((n, nmodes))
            hygro = jnp.zeros((n, nmodes))
            for l in range(vol.shape[1]):
                vaerosol = vaerosol + vol[:, l, :]
                hygro = hygro + vol[:, l, :] \
                    * W(spec_mask[l, :], spechygro[l, :], 0.0)[None, :]
            vok = vaerosol > 1.0e-30
            hygro = W(vok, hygro / W(vok, vaerosol, 1.0), 0.0)
            vaerosol = W(vok, vaerosol * rho[:, c, None], 0.0)
            naermod = 0.5 * (numg[:, c, :] + numg[:, cm, :]) \
                * rho[:, c, None]
            naermod = jnp.maximum(naermod, vaerosol * v2nhi[None, :])
            naermod = jnp.minimum(naermod, vaerosol * v2nlo[None, :])

            in_cloud = c < jb
            smax_f = jnp.zeros(n)
            mq = mact & (qcic[:, c] >= QSMALL)
            lamc_sf = W(mq & (lamc[:, c] != 0.0), lamc[:, c], 1.0)
            smax_f = W(in_cloud & mq,
                       ncic[:, c] / lamc_sf
                       * _gamma(2.0 + pgam[:, c]) / _gamma(1.0 + pgam[:, c]),
                       smax_f)
            mqr = mact & (qric[:, c] >= QSMALL)
            lamr_sf = W(mqr & (lamr[:, c] != 0.0), lamr[:, c], 1.0)
            smax_f = W(in_cloud & mqr,
                       smax_f + nric[:, c] / lamr_sf, smax_f)

            tair_s = W(mact, t[:, c], 250.0)
            rho_s = W(mact, rho[:, c], 1.0)
            wu_s = W(mact, wu[:, c], 1.0)
            fn = actdrop_mam_calc(wu_s, tair_s, rho_s, naermod,
                                  vaerosol, hygro, in_cloud,
                                  smax_f, ap)
            # sequential mode accumulation (Fortran loop order)
            nlsrc = jnp.zeros(n)
            for m in range(nmodes):
                nlsrc = nlsrc + fn[:, m] * naermod[:, m]
            dum2l = dum2l.at[:, c].set(W(mact, nlsrc, 0.0))

            # droplet activation rate (1997-2012)
            mrate = (act & (qcic[:, c] >= QSMALL) & (t[:, c] > 238.15)
                     & (c > jt + 2))
            npccn = W(mrate,
                      W(c == kqc, _divs(dum2l[:, c], deltat),
                        _divs(dum2l[:, c] - ncic[:, c], deltat)), 0.0)
            npccn = jnp.maximum(0.0, npccn)

            # ---- ice nucleation (2016-2152) ----
            es_c = svp_water(t[:, c])
            esi_c = svp_ice(t[:, c])
            qs_c = 0.622 * es_c / (ph[:, c] - (1.0 - 0.622) * es_c)
            qs_c = jnp.minimum(1.0, qs_c)
            qs_c = W(qs_c < 0.0, 1.0, qs_c)
            relhum_c = jnp.ones(n)

            mcold = act & (t[:, c] < TMELT)
            soot_num = 0.5 * (numg[:, cm, m_acc] + numg[:, c, m_acc]) \
                * rho[:, c] * 1.0e-6
            dmc = 0.5 * (mmrg[:, cm, l_dst, m_crs]
                         + mmrg[:, c, l_dst, m_crs])
            ssmc = 0.5 * (mmrg[:, cm, l_ncl, m_crs]
                          + mmrg[:, c, l_ncl, m_crs])
            so4mc = 0.5 * (mmrg[:, cm, l_so4, m_crs]
                           + mmrg[:, c, l_so4, m_crs])
            wght = dmc / W(dmc > 0.0, ssmc + dmc + so4mc, 1.0)
            dst_num = W(dmc > 0.0,
                        wght * (numg[:, cm, m_crs] + numg[:, c, m_crs])
                        * 0.5 * rho[:, c] * 1.0e-6, 0.0)
            dgnum_ait = 0.5 * (dgnumg[:, c, m_ait] + dgnumg[:, cm, m_ait])
            dga_pos = dgnum_ait > 0.0
            so4_num = W(dga_pos,
                        0.5 * (numg[:, cm, m_ait] + numg[:, c, m_ait])
                        * rho[:, c] * 1.0e-6
                        * (0.5 - 0.5 * _erf(_divs(
                            jnp.log(0.1e-6 / W(dga_pos, dgnum_ait, 1.0)),
                            2.0 ** 0.5 * np.log(sigmag_aitken)))),
                        0.0)
            so4_num = jnp.maximum(0.0, so4_num)
            soot_num_z = jnp.zeros(n)   # *** soot nucleation off ***

            mnuc = mcold & (wu[:, c] < 4.0)
            nuci = nucleati_conv(
                W(mnuc, wu[:, c], 0.5), W(mnuc, t[:, c], 250.0),
                relhum_c, jnp.ones(n), W(mnuc, qcic[:, c], 0.0),
                W(mnuc, rho[:, c], 1.0), W(mnuc, so4_num, 0.0),
                W(mnuc, dst_num, 0.0), soot_num_z)
            dum2i = dum2i.at[:, c].set(W(mcold & mnuc, nuci,
                                         W(mcold, dum2i[:, c], 0.0)))

            mice = (act & (dum2i[:, c] > 0.0) & (t[:, c] < TMELT)
                    & (relhum_c * es_c / esi_c > 1.05) & (c > jt + 1))
            nnuccd = W(mice,
                       W(c == kqi, _divs(dum2i[:, c], deltat),
                         _divs(dum2i[:, c] - niic[:, c], deltat)), 0.0)
            nnuccd = jnp.maximum(nnuccd, 0.0)
            mnuccd = nnuccd * mi0

            # ---- Bergeron process, Rotstayn (2169-2191) ----
            mbg = (act & (t[:, c] <= 273.15) & (t[:, c] > 233.15)
                   & (qiic[:, c] > 0.5e-6) & (qcic[:, c] > QSMALL))
            a_pr = LS_B * (LS_B / (RV_B * t[:, c]) - 1.0) \
                / (KA_B * t[:, c])
            b_pr = RV_B * t[:, c] * ph[:, c] / (2.21 * esi_c)
            cpvd = 65.2 * W(mbg, niic[:, c], 0.0) ** 0.5 \
                * (es_c - esi_c) / ((a_pr + b_pr) * esi_c)
            dqi = jnp.maximum(
                0.0, (0.5 * cpvd * deltat
                      + W(mbg, qiic[:, c], 0.0) ** 0.5) ** 2.0
                - qiic[:, c])
            prb = W(mbg, _divs(jnp.minimum(qcic[:, c], dqi), deltat), 0.0)
            # prb/(qcic/ncic): IEEE x/0=inf, y/inf=0 when ncic=0
            nprb = W(mbg, prb / (W(mbg, qcic[:, c], 1.0)
                                 / W(mbg, ncic[:, c], 1.0)), 0.0)

            # ---- het. freezing of cloud water (2196-2301) ----
            mhf = (act & (qcic[:, c] >= QSMALL) & (ncic[:, c] > 0.0)
                   & (ntaerh_c > 0.0) & (t[:, c] <= 268.15)
                   & (t[:, c] > 238.15))
            lamc_s5 = W(mhf & (lamc[:, c] != 0.0), lamc[:, c], 1.0)
            g_7pg = _gamma(7.0 + pgam[:, c])
            g_pg4b = _gamma(pgam[:, c] + 4.0)
            mnuccc = W(mhf, PI * PI / 36.0 * RHOW * cdist1[:, c]
                       * g_7pg * BIMM
                       * (jnp.exp(AIMM * (273.15 - t[:, c])) - 1.0)
                       / lamc_s5 ** 3.0 / lamc_s5 ** 3.0, 0.0)
            nnuccc = W(mhf, PI / 6.0 * cdist1[:, c] * g_pg4b * BIMM
                       * (jnp.exp(AIMM * (273.15 - t[:, c])) - 1.0)
                       / lamc_s5 ** 3.0, 0.0)

            tcnt = (270.16 - t[:, c]) ** 1.3
            viscosity = 1.8e-5 * _divs(t[:, c], 298.0) ** 0.85
            mfp = 2.0 * viscosity / (ph[:, c] * jnp.sqrt(
                8.0 * 28.96e-3 / (PI * 8.314409 * t[:, c])))
            slip1 = 1.0 + _divs(mfp, RN_DST1) * (
                1.257 + 0.4 * jnp.exp(-(1.1 * RN_DST1 / mfp)))
            slip2 = 1.0 + _divs(mfp, RN_DST2) * (
                1.257 + 0.4 * jnp.exp(-(1.1 * RN_DST2 / mfp)))
            slip3 = 1.0 + _divs(mfp, RN_DST3) * (
                1.257 + 0.4 * jnp.exp(-(1.1 * RN_DST3 / mfp)))
            slip4 = 1.0 + _divs(mfp, RN_DST4) * (
                1.257 + 0.4 * jnp.exp(-(1.1 * RN_DST4 / mfp)))
            dfaer1 = 1.381e-23 * t[:, c] * slip1 \
                / (6.0 * PI * viscosity * RN_DST1)
            dfaer2 = 1.381e-23 * t[:, c] * slip2 \
                / (6.0 * PI * viscosity * RN_DST2)
            dfaer3 = 1.381e-23 * t[:, c] * slip3 \
                / (6.0 * PI * viscosity * RN_DST3)
            dfaer4 = 1.381e-23 * t[:, c] * slip4 \
                / (6.0 * PI * viscosity * RN_DST4)
            nacon3 = W(dmc > 0.0, dst_num * tcnt * 1.0e6, 0.0)
            mnucct = W(mhf, (dfaer3 * nacon3) * PI * PI / 3.0 * RHOW
                       * cdist1[:, c] * _gamma(pgam[:, c] + 5.0)
                       / lamc_s5 ** 4.0, 0.0)
            nnucct = W(mhf, (dfaer3 * nacon3) * 2.0 * PI
                       * cdist1[:, c] * _gamma(pgam[:, c] + 2.0)
                       / lamc_s5, 0.0)
            _ = (dfaer1, dfaer2, dfaer4)  # nacon1/2/4 = 0 (modal)

            # ---- homogeneous freezing of cloud water (2304-2319) ----
            mhom = act & (t[:, c] < 233.15) & (qc[:, c] > 0.0)
            dum = xlf / cp * qc[:, c]
            partial = mhom & (t[:, c] + dum > 233.15)
            dum = W(partial,
                    jnp.clip(_divs(-(t[:, c] - 233.15) * cp, xlf)
                             / W(mhom, qc[:, c], 1.0), 0.0, 1.0),
                    1.0)
            fholm = fholm.at[:, c].set(
                W(mhom, muc * dum * qc[:, c], fholm[:, c]))
            fholn = fholn.at[:, c].set(
                W(mhom, muc * dum * nc[:, c], fholn[:, c]))

            # =======================================================
            # conservation checks (2340-2552)
            # =======================================================
            arc = arcf[:, c]
            qce = muc * qc[:, c] - fholm[:, c] + dzc * cmel[:, cm]
            dum = arc * (pra + prc + prb + mnuccc + mnucct + msacwi
                         + psacws + psacwg + pgsacw + qmultg
                         + psacwi) * dzc
            neg = act & (qce < 0.0)
            over = act & ~neg & (dum > qce)
            ratio = qce / W(over, dum, 1.0) * omsm
            prc = W(neg, 0.0, prc)
            pra = W(neg, 0.0, pra)
            prb = W(neg, 0.0, prb)
            mnuccc = W(neg, 0.0, mnuccc)
            mnucct = W(neg, 0.0, mnucct)
            msacwi = W(neg, 0.0, msacwi)
            psacws = W(neg, 0.0, psacws)
            psacwg = W(neg, 0.0, psacwg)
            pgsacw = W(neg, 0.0, pgsacw)
            qmultg = W(neg, 0.0, qmultg)
            psacwi = W(neg, 0.0, psacwi)
            prc = W(over, prc * ratio, prc)
            pra = W(over, pra * ratio, pra)
            prb = W(over, prb * ratio, prb)
            mnuccc = W(over, mnuccc * ratio, mnuccc)
            mnucct = W(over, mnucct * ratio, mnucct)
            msacwi = W(over, msacwi * ratio, msacwi)
            psacws = W(over, psacws * ratio, psacws)
            psacwg = W(over, psacwg * ratio, psacwg)
            pgsacw = W(over, pgsacw * ratio, pgsacw)
            qmultg = W(over, qmultg * ratio, qmultg)
            psacwi = W(over, psacwi * ratio, psacwi)

            nce = (muc * nc[:, c] - fholn[:, c]
                   + (arc * npccn * mtimec) * dzc)
            dum = arc * dzc * (nprc1 + npra + nnuccc + nnucct
                               + npsacws + nprb + npsacwg + npsacwi)
            neg = act & (nce < 0.0)
            over = act & ~neg & (dum > nce)
            ratio = nce / W(over, dum, 1.0) * omsm
            nprc1 = W(neg, 0.0, W(over, nprc1 * ratio, nprc1))
            npra = W(neg, 0.0, W(over, npra * ratio, npra))
            nnuccc = W(neg, 0.0, W(over, nnuccc * ratio, nnuccc))
            nnucct = W(neg, 0.0, W(over, nnucct * ratio, nnucct))
            npsacws = W(neg, 0.0, W(over, npsacws * ratio, npsacws))
            nprb = W(neg, 0.0, W(over, nprb * ratio, nprb))
            npsacwg = W(neg, 0.0, W(over, npsacwg * ratio, npsacwg))
            npsacwi = W(neg, 0.0, W(over, npsacwi * ratio, npsacwi))

            qre = muc * qr[:, c] + dzc * (pra + prc) * arc
            dum = arc * dzc * (pracs + mnuccr - prf[:, c] + pracg
                               + pgracs + piacr + piacrs + qmultrg)
            neg = act & (qre < 0.0)
            over = act & ~neg & (dum > qre)
            ratio = qre / W(over, dum, 1.0) * omsm
            prf = prf.at[:, c].set(
                W(neg, 0.0, W(over, prf[:, c] * ratio, prf[:, c])))
            pracs = W(neg, 0.0, W(over, pracs * ratio, pracs))
            mnuccr = W(neg, 0.0, W(over, mnuccr * ratio, mnuccr))
            pracg = W(neg, 0.0, W(over, pracg * ratio, pracg))
            pgracs = W(neg, 0.0, W(over, pgracs * ratio, pgracs))
            piacr = W(neg, 0.0, W(over, piacr * ratio, piacr))
            piacrs = W(neg, 0.0, W(over, piacrs * ratio, piacrs))
            qmultrg = W(neg, 0.0, W(over, qmultrg * ratio, qmultrg))

            # top-of-plume qr fix (2424-2434)
            mtop_r = act & (cm == jt + 1)
            dum = (pra + prc - pracs - mnuccr - pracg - pgracs - piacr
                   - piacrs - qmultrg)
            fix = mtop_r & (dum < 0.0)
            mn_pos = fix & (mnuccr > 0.0)
            mnuccr = W(mn_pos, jnp.maximum(mnuccr + _divs(dum, omsm), 0.0),
                       mnuccr)
            pracg = W(fix & ~mn_pos,
                      jnp.maximum(pracg + _divs(dum, omsm), 0.0), pracg)

            nre = muc * nr[:, c] + nprc * arc * dzc
            dum = arc * dzc * (npracs + nnuccr - nragg - pnrf[:, c]
                               + niacr + niacrs + npracg - ngracs)
            neg = act & (nre < 0.0)
            over = act & ~neg & (dum > nre)
            ratio = nre / W(over, dum, 1.0) * omsm
            npracs = W(neg, 0.0, W(over, npracs * ratio, npracs))
            nnuccr = W(neg, 0.0, W(over, nnuccr * ratio, nnuccr))
            nragg = W(neg, 0.0, W(over, nragg * ratio, nragg))
            pnrf = pnrf.at[:, c].set(
                W(neg, 0.0, W(over, pnrf[:, c] * ratio, pnrf[:, c])))
            niacr = W(neg, 0.0, W(over, niacr * ratio, niacr))
            niacrs = W(neg, 0.0, W(over, niacrs * ratio, niacrs))
            npracg = W(neg, 0.0, W(over, npracg * ratio, npracg))
            ngracs = W(neg, 0.0, W(over, ngracs * ratio, ngracs))

            qie = (muc * qi[:, c] + fholm[:, c]
                   + dzc * (cmei[:, cm]
                            + (mnuccc + mnucct + msacwi + prb + qmultg
                               + qmultrg + psacwi) * arc))
            dum = arc * (prci + prai + praci + pracis) * dzc
            neg = act & (qie < 0.0)
            over = act & ~neg & (dum > qie)
            ratio = qie / W(over, dum, 1.0) * omsm
            prci = W(neg, 0.0, W(over, prci * ratio, prci))
            prai = W(neg, 0.0, W(over, prai * ratio, prai))
            praci = W(neg, 0.0, W(over, praci * ratio, praci))
            pracis = W(neg, 0.0, W(over, pracis * ratio, pracis))

            nie = (muc * ni[:, c] + fholn[:, c]
                   + dzc * (nnuccd * mtime * arc
                            + (nnuccc + nnucct + nmultg + nmultrg)
                            * arc))
            dum = arc * dzc * (-nsacwi + nprci + nprai + niacr + niacrs)
            neg = act & (nie < 0.0)
            over = act & ~neg & (dum > nie)
            ratio = nie / W(over, dum, 1.0) * omsm
            nsacwi = W(neg, 0.0, W(over, nsacwi * ratio, nsacwi))
            nprci = W(neg, 0.0, W(over, nprci * ratio, nprci))
            nprai = W(neg, 0.0, W(over, nprai * ratio, nprai))
            niacr = W(neg, 0.0, W(over, niacr * ratio, niacr))
            niacrs = W(neg, 0.0, W(over, niacrs * ratio, niacrs))

            qnie = (muc * qni[:, c]
                    + dzc * ((prai + psacws + prci + pracs + piacrs
                              + pracis) * arc))
            dum = arc * dzc * (-psf[:, c] + psacr)
            neg = act & (qnie < 0.0)
            over = act & ~neg & (dum > qnie)
            ratio = qnie / W(over, dum, 1.0) * omsm
            psf = psf.at[:, c].set(
                W(neg, 0.0, W(over, psf[:, c] * ratio, psf[:, c])))
            psacr = W(neg, 0.0, W(over, psacr * ratio, psacr))

            nse = muc * ns[:, c] + dzc * (nprci + niacrs) * arc
            dum = arc * dzc * (-nsagg - pnsf[:, c] + nscng + ngracs)
            neg = act & (nse < 0.0)
            over = act & ~neg & (dum > nse)
            ratio = nse / W(over, dum, 1.0) * omsm
            nsagg = W(neg, 0.0, W(over, nsagg * ratio, nsagg))
            pnsf = pnsf.at[:, c].set(
                W(neg, 0.0, W(over, pnsf[:, c] * ratio, pnsf[:, c])))
            nscng = W(neg, 0.0, W(over, nscng * ratio, nscng))
            ngracs = W(neg, 0.0, W(over, ngracs * ratio, ngracs))

            qge = (muc * qg[:, c]
                   + dzc * (pracg + psacwg + pgsacw + pgracs + mnuccr
                            + piacr + praci + psacr) * arc)
            dum = arc * dzc * (-pgf[:, c])
            neg = act & (qge < 0.0)
            over = act & ~neg & (dum > qge)
            ratio = qge / W(over, dum, 1.0) * omsm
            pgf = pgf.at[:, c].set(
                W(neg, 0.0, W(over, pgf[:, c] * ratio, pgf[:, c])))

            nge = (muc * ng[:, c]
                   + dzc * (nscng + ngracs + nnuccr + niacr) * arc)
            dum = arc * dzc * (-pngf[:, c])
            neg = act & (nge < 0.0)
            over = act & ~neg & (dum > nge)
            ratio = nge / W(over, dum, 1.0) * omsm
            pngf = pngf.at[:, c].set(
                W(neg, 0.0, W(over, pngf[:, c] * ratio, pngf[:, c])))

            # =======================================================
            # tendencies (2559-2671)
            # =======================================================
            mten = act & (c <= kqc)
            qctend = qctend.at[:, c].set(W(mten,
                (-pra - prc - prb - mnuccc - mnucct - msacwi - psacws)
                - psacwg - pgsacw - qmultg - psacwi, W(act, 0.0,
                                                       qctend[:, c])))
            qitend = qitend.at[:, c].set(W(mten,
                (prb + mnuccc + mnucct + msacwi - prci - prai) - praci
                - pracis + qmultg + qmultrg + psacwi,
                W(act, 0.0, qitend[:, c])))
            qrtend_c = (pra + prc - pracs - mnuccr - pracg - pgracs
                        - piacr - piacrs - qmultrg)
            qnitend_c = ((prai + psacws + prci) + pracs - psacr
                         + piacrs + pracis)
            qgtend_c = (pracg + psacwg + pgsacw + pgracs + mnuccr
                        + piacr + praci + psacr)
            neg_sum = mten & ((qnitend_c + qrtend_c + qgtend_c) < 0.0)
            dum = _divs(qnitend_c + qrtend_c + qgtend_c, omsm)
            qrtend_c = W(neg_sum, qrtend_c - dum, qrtend_c)
            qmultrg = W(neg_sum, qmultrg + dum, qmultrg)
            qitend = qitend.at[:, c].set(
                W(neg_sum, qitend[:, c] + dum, qitend[:, c]))
            qrtend = qrtend.at[:, c].set(W(mten, qrtend_c,
                                           W(act, 0.0, qrtend[:, c])))
            qnitend = qnitend.at[:, c].set(
                W(mten, qnitend_c, W(act, 0.0, qnitend[:, c])))
            qgtend = qgtend.at[:, c].set(W(mten, qgtend_c,
                                           W(act, 0.0, qgtend[:, c])))

            nctend = nctend.at[:, c].set(W(mten,
                npccn * mtimec + (-nnuccc - nnucct - npsacws - npra
                                  - nprc1 - nprb - npsacwg - npsacwi),
                W(act, 0.0, nctend[:, c])))
            nitend = nitend.at[:, c].set(W(mten,
                nnuccd * mtime + (nnuccc + nnucct + nsacwi - nprci
                                  - nprai) - niacr - niacrs + nmultg
                + nmultrg, W(act, 0.0, nitend[:, c])))
            nstend = nstend.at[:, c].set(W(mten,
                nsagg + nprci - nscng - ngracs + niacrs,
                W(act, 0.0, nstend[:, c])))
            nrtend = nrtend.at[:, c].set(W(mten,
                nprc + (-npracs - nnuccr + nragg) - niacr - niacrs
                - npracg - ngracs, W(act, 0.0, nrtend[:, c])))
            ngtend = ngtend.at[:, c].set(W(mten,
                nscng + ngracs + nnuccr + niacr,
                W(act, 0.0, ngtend[:, c])))

            # diagnostic outputs at k-1 (2600-2659), mten lanes only
            def dset(nm, val):
                diag[nm] = diag[nm].at[:, cm].set(
                    W(mten, val, diag[nm][:, cm]))
            dset("autolm", -prc * arc)
            dset("accrlm", -pra * arc)
            dset("bergnm", -prb * arc)
            dset("fhtimm", -mnuccc * arc)
            dset("fhtctm", -mnucct * arc)
            dset("hmpim", -msacwi * arc)
            dset("accslm", -psacws * arc)
            dset("fhmlm", -fholm[:, c] / dzc)
            dset("autoln", -nprc1 * arc)
            dset("accrln", -npra * arc)
            dset("bergnn", -nprb * arc)
            dset("fhtimn", -nnuccc * arc)
            dset("fhtctn", -nnucct * arc)
            dset("accsln", -npsacws * arc)
            dset("activn", npccn * mtimec * arc)
            dset("fhmln", -fholn[:, c] / dzc)
            dset("autoim", -prci * arc)
            dset("accsim", -prai * arc)
            dset("nuclin", nnuccd * mtime * arc)
            dset("autoin", -nprci * arc)
            dset("accsin", -nprai * arc)
            dset("hmpin", nsacwi * arc)
            dset("accgrm", -pracg * arc)
            dset("accglm", -psacwg * arc)
            dset("accgslm", -pgsacw * arc)
            dset("accgsrm", -pgracs * arc)
            dset("accgirm", -piacr * arc)
            dset("accgrim", -praci * arc)
            dset("accgrsm", -psacr * arc)
            dset("accgsln", -nscng * arc)
            dset("accgsrn", -ngracs * arc)
            dset("accgirn", -niacr * arc)
            dset("accsrim", -pracis * arc)
            dset("acciglm", -qmultg * arc)
            dset("accigrm", qmultrg * arc)
            dset("accsirm", -piacrs * arc)
            dset("accigln", nmultg * arc)
            dset("accigrn", nmultrg * arc)
            dset("accsirn", -niacrs * arc)
            dset("accgln", -npsacwg * arc)
            dset("accgrn", -npracg * arc)
            dset("accilm", -psacwi * arc)
            dset("acciln", -npsacwi * arc)
            dset("fallrm", prf[:, c] * arc)
            dset("fallsm", psf[:, c] * arc)
            dset("fallgm", pgf[:, c] * arc)
            dset("fallrn", pnrf[:, c] * arc)
            dset("fallsn", pnsf[:, c] * arc)
            dset("fallgn", pngf[:, c] * arc)

            # =======================================================
            # vertical integration (2677-2861)
            # =======================================================
            mig = act & (c <= kqi)
            mucm_s = W(act, mucm, 1.0)
            qg_new = 1.0 / mucm_s * (muc * qg[:, c]
                                     + dzc * (qgtend[:, c] + pgf[:, c])
                                     * arc)
            ng_new = 1.0 / mucm_s * (muc * ng[:, c]
                                     + dzc * (ngtend[:, c] + pngf[:, c])
                                     * arc)
            qg = qg.at[:, cm].set(W(mig, qg_new, W(act, 0.0, qg[:, cm])))
            ng = ng.at[:, cm].set(
                W(mig, jnp.maximum(ng_new, 0.0), W(act, 0.0, ng[:, cm])))
            gz = act & (qg[:, cm] <= 0.0)
            qg = qg.at[:, cm].set(W(gz, 0.0, qg[:, cm]))
            ng = ng.at[:, cm].set(W(gz, 0.0, ng[:, cm]))

            den_s = W(act, mucm + dzc * du[:, cm], 1.0)
            qni_new = 1.0 / den_s * (muc * qni[:, c]
                                     + dzc * (qnitend[:, c] + psf[:, c])
                                     * arc)
            ns_new = 1.0 / den_s * (muc * ns[:, c]
                                    + dzc * (nstend[:, c] + pnsf[:, c])
                                    * arc)
            qni = qni.at[:, cm].set(W(mig, qni_new,
                                      W(act, 0.0, qni[:, cm])))
            ns = ns.at[:, cm].set(
                W(mig, jnp.maximum(ns_new, 0.0), W(act, 0.0, ns[:, cm])))
            qnide = qnide.at[:, c].set(W(mig, qni[:, cm], qnide[:, c]))
            nsde = nsde.at[:, c].set(W(mig, ns[:, cm], nsde[:, c]))
            diag["dsfm"] = diag["dsfm"].at[:, cm].set(
                W(act, -du[:, cm] * qnide[:, c], diag["dsfm"][:, cm]))
            diag["dsfn"] = diag["dsfn"].at[:, cm].set(
                W(act, -du[:, cm] * nsde[:, c], diag["dsfn"][:, cm]))
            sz = act & (qni[:, cm] <= 0.0)
            qni = qni.at[:, cm].set(W(sz, 0.0, qni[:, cm]))
            ns = ns.at[:, cm].set(W(sz, 0.0, ns[:, cm]))

            mrc = act & (c <= kqc)
            qr_new = 1.0 / mucm_s * (muc * qr[:, c]
                                     + dzc * (qrtend[:, c] + prf[:, c])
                                     * arc)
            nr_new = 1.0 / mucm_s * (muc * nr[:, c]
                                     + dzc * (nrtend[:, c] + pnrf[:, c])
                                     * arc)
            qr = qr.at[:, cm].set(W(mrc, qr_new, W(act, 0.0, qr[:, cm])))
            nr = nr.at[:, cm].set(
                W(mrc, jnp.maximum(nr_new, 0.0), W(act, 0.0, nr[:, cm])))
            rz = act & (qr[:, cm] <= 0.0)
            qr = qr.at[:, cm].set(W(rz, 0.0, qr[:, cm]))
            nr = nr.at[:, cm].set(W(rz, 0.0, nr[:, cm]))

            # homogeneous rain freezing at k-1 (2738-2768)
            mfrz = act & (t[:, cm] < 233.15) & (qr[:, cm] > 0.0)
            dum = xlf / cp * qr[:, cm]
            partial = mfrz & (t[:, cm] + dum > 233.15)
            dum = W(partial,
                    jnp.clip(_divs(-(t[:, cm] - 233.15) * cp, xlf)
                             / W(mfrz, qr[:, cm], 1.0), 0.0, 1.0),
                    1.0)
            fhmrm_new = -mucm * dum * qr[:, cm] / dzc
            fhmrm = fhmrm.at[:, cm].set(W(mfrz, fhmrm_new,
                                          fhmrm[:, cm]))
            mtfix = (mfrz & (cm == jt + 1)
                     & ((qrtend[:, c] * arc + fhmrm[:, cm]) < 0.0))
            fhmrm = fhmrm.at[:, cm].set(
                W(mtfix, -qrtend[:, c] * arc, fhmrm[:, cm]))
            dum = W(mtfix,
                    -fhmrm[:, cm] * dzc / mucm_s
                    / W(mfrz, qr[:, cm], 1.0), dum)
            qg = qg.at[:, cm].set(
                W(mfrz, qg[:, cm] + dum * qr[:, cm], qg[:, cm]))
            ng = ng.at[:, cm].set(
                W(mfrz, jnp.maximum(ng[:, cm] + dum * nr[:, cm], 0.0),
                  ng[:, cm]))
            qr = qr.at[:, cm].set(
                W(mfrz, (1.0 - dum) * qr[:, cm], qr[:, cm]))
            nr = nr.at[:, cm].set(
                W(mfrz, jnp.maximum((1.0 - dum) * nr[:, cm], 0.0),
                  nr[:, cm]))

            # snow detrainment adjustment (2771-2782)
            dum = jnp.minimum(
                (qnitend[:, c] + qrtend[:, c] + qgtend[:, c]) * arc,
                (qnitend[:, c] + qgtend[:, c]) * arc - fhmrm[:, cm])
            madj = (act & ((dum + diag["dsfm"][:, cm]) < 0.0)
                    & (diag["dsfm"][:, cm] != 0.0))
            dsfm_old = W(madj, diag["dsfm"][:, cm], 1.0)
            diag["dsfn"] = diag["dsfn"].at[:, cm].set(
                W(madj, diag["dsfn"][:, cm] * (-dum * omsm) / dsfm_old,
                  diag["dsfn"][:, cm]))
            diag["dsfm"] = diag["dsfm"].at[:, cm].set(
                W(madj, -dum * omsm, diag["dsfm"][:, cm]))
            du_s = W(act & (du[:, cm] != 0.0), du[:, cm],
                     W(madj, jnp.inf, 1.0))
            qnide = qnide.at[:, c].set(
                W(madj, -diag["dsfm"][:, cm] / du_s, qnide[:, c]))
            nsde = nsde.at[:, c].set(
                W(madj, -diag["dsfn"][:, cm] / du_s, nsde[:, c]))
            qni = qni.at[:, cm].set(
                W(madj, qni[:, cm] + dzc * du[:, cm]
                  * (qni[:, cm] - qnide[:, c]) / mucm_s, qni[:, cm]))
            ns = ns.at[:, cm].set(
                W(madj, jnp.maximum(
                    ns[:, cm] + dzc * du[:, cm]
                    * (ns[:, cm] - nsde[:, c]) / mucm_s, 0.0),
                  ns[:, cm]))

            rprd = rprd.at[:, cm].set(W(act,
                (qnitend[:, c] + qrtend[:, c] + qgtend[:, c]) * arc
                + diag["dsfm"][:, cm], rprd[:, cm]))
            sprd = sprd.at[:, cm].set(W(act,
                (qnitend[:, c] + qgtend[:, c]) * arc - fhmrm[:, cm]
                + diag["dsfm"][:, cm], sprd[:, cm]))

            # cloud water integration (2788-2817)
            qc_new = (muc * qc[:, c] - fholm[:, c]
                      + dzc * qctend[:, c] * arc
                      + dzc * cmel[:, cm]) / den_s
            nc_new = (muc * nc[:, c] - fholn[:, c]
                      + dzc * nctend[:, c] * arc) / den_s
            qc = qc.at[:, cm].set(W(mrc, qc_new, W(act, 0.0, qc[:, cm])))
            qcde = qcde.at[:, c].set(W(mrc, qc[:, cm], qcde[:, c]))
            nc = nc.at[:, cm].set(
                W(mrc, jnp.maximum(nc_new, 0.0), W(act, 0.0, nc[:, cm])))
            ncde = ncde.at[:, c].set(W(mrc, nc[:, cm], ncde[:, c]))
            diag["dlfm"] = diag["dlfm"].at[:, cm].set(
                W(act, -du[:, cm] * qcde[:, c], diag["dlfm"][:, cm]))
            diag["dlfn"] = diag["dlfn"].at[:, cm].set(
                W(act, -du[:, cm] * ncde[:, c], diag["dlfn"][:, cm]))
            cz = act & (qc[:, cm] <= 0.0)
            qc = qc.at[:, cm].set(W(cz, 0.0, qc[:, cm]))
            nc = nc.at[:, cm].set(W(cz, 0.0, nc[:, cm]))

            # cloud ice integration (2820-2848)
            qi_new = (muc * qi[:, c] + fholm[:, c]
                      + dzc * qitend[:, c] * arc
                      + dzc * cmei[:, cm]) / den_s
            ni_new = (muc * ni[:, c] + fholn[:, c]
                      + dzc * nitend[:, c] * arc) / den_s
            qi = qi.at[:, cm].set(W(mig, qi_new, W(act, 0.0, qi[:, cm])))
            qide = qide.at[:, c].set(W(mig, qi[:, cm], qide[:, c]))
            ni = ni.at[:, cm].set(
                W(mig, jnp.maximum(ni_new, 0.0), W(act, 0.0, ni[:, cm])))
            nide = nide.at[:, c].set(W(mig, ni[:, cm], nide[:, c]))
            diag["difm"] = diag["difm"].at[:, cm].set(
                W(act, -du[:, cm] * qide[:, c], diag["difm"][:, cm]))
            diag["difn"] = diag["difn"].at[:, cm].set(
                W(act, -du[:, cm] * nide[:, c], diag["difn"][:, cm]))
            iz = act & (qi[:, cm] <= 0.0)
            qi = qi.at[:, cm].set(W(iz, 0.0, qi[:, cm]))
            ni = ni.at[:, cm].set(W(iz, 0.0, ni[:, cm]))

            frz = frz.at[:, cm].set(W(act,
                cmei[:, cm] + arc * (prb + mnuccc + mnucct + msacwi
                                     + pracs + mnuccr + psacws + pracg
                                     + psacwg + pgsacw + pgracs + piacr
                                     + piacrs + qmultg + qmultrg
                                     + psacwi)
                - diag["fhmlm"][:, cm] - fhmrm[:, cm], frz[:, cm]))

            # ---- post-integration size distributions at k-1 ----
            # cloud ice (2877-2939) with niadj redistribution
            niorg = ni[:, cm]
            mi2 = act & (qi[:, cm] >= QSMALL)
            ni = ni.at[:, cm].set(
                W(mi2, jnp.minimum(ni[:, cm], qi[:, cm] * 1.0e20),
                  ni[:, cm]))
            qi_s = W(mi2, qi[:, cm], 1.0)
            lami_m = (mp["g_1di"] * ci * W(mi2, ni[:, cm], 1.0)
                      / qi_s) ** (1.0 / di)
            lo = mi2 & (lami_m < lammin_i)
            hi = mi2 & (lami_m > lammax_i)
            lami_m = W(lo, lammin_i, W(hi, lammax_i, lami_m))
            n0i_m = W(lo | hi,
                      _divs(lami_m ** (di + 1.0) * qi_s,
                            ci * mp["g_1di"]),
                      W(mi2, ni[:, cm], 0.0) * lami_m)
            ni = ni.at[:, cm].set(W(lo | hi, n0i_m / lami_m, ni[:, cm]))
            lami = lami.at[:, cm].set(W(mi2, lami_m,
                                        W(act, 0.0, lami[:, cm])))
            n0i = n0i.at[:, cm].set(W(mi2, n0i_m,
                                      W(act, 0.0, n0i[:, cm])))
            nide = nide.at[:, c].set(W(act, ni[:, cm], nide[:, c]))
            diag["difn"] = diag["difn"].at[:, cm].set(
                W(act, -du[:, cm] * nide[:, c], diag["difn"][:, cm]))
            niadj_c = (ni[:, cm] - niorg) * mucm / dzc
            m_neg = act & (niadj_c < 0.0)
            total = (diag["nuclin"][:, cm] - diag["fhtimn"][:, cm]
                     - diag["fhtctn"][:, cm] - diag["fhmln"][:, cm]
                     + diag["hmpin"][:, cm])
            tz = total != 0.0
            for nm in ("nuclin", "fhtimn", "fhtctn", "fhmln", "hmpin"):
                v = diag[nm][:, cm]
                v_new = W(tz, v + v * niadj_c / W(tz, total, 1.0),
                          v + _divs(niadj_c, 5.0))
                diag[nm] = diag[nm].at[:, cm].set(W(m_neg, v_new, v))
            m_pos = act & (niadj_c > 0.0)
            total = diag["autoin"][:, cm] + diag["accsin"][:, cm]
            tz = total != 0.0
            for nm in ("autoin", "accsin"):
                v = diag[nm][:, cm]
                v_new = W(tz, v + v * niadj_c / W(tz, total, 1.0),
                          v + niadj_c / 2.0)
                diag[nm] = diag[nm].at[:, cm].set(W(m_pos, v_new, v))

            # cloud water (2943-3006) with ncadj redistribution
            ncorg = nc[:, cm]
            mc2 = act & (qc[:, cm] >= QSMALL)
            nc = nc.at[:, cm].set(
                W(mc2, jnp.minimum(nc[:, cm], qc[:, cm] * 1.0e20),
                  nc[:, cm]))
            pg_m = 0.0005714 * (_divs(W(mc2, nc[:, cm], 0.0), 1.0e6)
                                / rho[:, cm]) + 0.2714
            pg_m = 1.0 / (pg_m ** 2.0) - 1.0
            pg_m = jnp.clip(pg_m, 2.0, 15.0)
            pgam = pgam.at[:, cm].set(W(mc2, pg_m, pgam[:, cm]))
            qc_s = W(mc2, qc[:, cm], 1.0)
            g_pg4m = _gamma(pg_m + 4.0)
            g_pg1m = _gamma(pg_m + 1.0)
            lamc_m = (PI / 6.0 * RHOW * W(mc2, nc[:, cm], 1.0)
                      * g_pg4m / (qc_s * g_pg1m)) ** (1.0 / 3.0)
            lammin_m = _divs(pg_m + 1.0, 40.0e-6)
            lammax_m = _divs(pg_m + 1.0, 1.0e-6)
            lo = mc2 & (lamc_m < lammin_m)
            hi = mc2 & (lamc_m > lammax_m)
            lamc_m = W(lo, lammin_m, W(hi, lammax_m, lamc_m))
            nc = nc.at[:, cm].set(
                W(lo | hi, 6.0 * lamc_m ** 3.0 * qc_s * g_pg1m
                  / (PI * RHOW * g_pg4m), nc[:, cm]))
            lamc = lamc.at[:, cm].set(W(mc2, lamc_m,
                                        W(act, 0.0, lamc[:, cm])))
            cdist1 = cdist1.at[:, cm].set(
                W(mc2, nc[:, cm] / g_pg1m, W(act, 0.0, cdist1[:, cm])))
            ncde = ncde.at[:, c].set(W(act, nc[:, cm], ncde[:, c]))
            diag["dlfn"] = diag["dlfn"].at[:, cm].set(
                W(act, -du[:, cm] * ncde[:, c], diag["dlfn"][:, cm]))
            ncadj_c = (nc[:, cm] - ncorg) * mucm / dzc
            m_neg = act & (ncadj_c < 0.0)
            diag["activn"] = diag["activn"].at[:, cm].set(
                W(m_neg, diag["activn"][:, cm] + ncadj_c,
                  diag["activn"][:, cm]))
            m_pos = act & (ncadj_c > 0.0)
            total = (diag["autoln"][:, cm] + diag["accrln"][:, cm]
                     + diag["bergnn"][:, cm] + diag["accsln"][:, cm])
            tz = total != 0.0
            for nm in ("autoln", "accrln", "bergnn", "accsln"):
                v = diag[nm][:, cm]
                v_new = W(tz, v + v * ncadj_c / W(tz, total, 1.0),
                          v + ncadj_c / 4.0)
                diag[nm] = diag[nm].at[:, cm].set(W(m_pos, v_new, v))

            # transport terms (3008-3026)
            diag["trspcm"] = diag["trspcm"].at[:, cm].set(W(act,
                (muc * qc[:, c] - mucm * qc[:, cm]) / dzc,
                diag["trspcm"][:, cm]))
            diag["trspcn"] = diag["trspcn"].at[:, cm].set(W(act,
                (muc * nc[:, c] - mucm * nc[:, cm]) / dzc,
                diag["trspcn"][:, cm]))
            diag["trspim"] = diag["trspim"].at[:, cm].set(W(act,
                (muc * qi[:, c] - mucm * qi[:, cm]) / dzc,
                diag["trspim"][:, cm]))
            diag["trspin"] = diag["trspin"].at[:, cm].set(W(act,
                (muc * ni[:, c] - mucm * ni[:, cm]) / dzc,
                diag["trspin"][:, cm]))
            mtt = act & (cm == jt + 1)
            if bool(jnp.any(mtt)) and c - 2 >= 0:
                cm2 = c - 2
                diag["trspcm"] = diag["trspcm"].at[:, cm2].set(
                    W(mtt, mucm * qc[:, cm] / dz[:, cm],
                      diag["trspcm"][:, cm2]))
                diag["trspcn"] = diag["trspcn"].at[:, cm2].set(
                    W(mtt, mucm * nc[:, cm] / dz[:, cm],
                      diag["trspcn"][:, cm2]))
                diag["trspim"] = diag["trspim"].at[:, cm2].set(
                    W(mtt, mucm * qi[:, cm] / dz[:, cm],
                      diag["trspim"][:, cm2]))
                diag["trspin"] = diag["trspin"].at[:, cm2].set(
                    W(mtt, mucm * ni[:, cm] / dz[:, cm],
                      diag["trspin"][:, cm2]))
                qcde = qcde.at[:, cm].set(W(mtt, qc[:, cm], qcde[:, cm]))
                ncde = ncde.at[:, cm].set(W(mtt, nc[:, cm], ncde[:, cm]))
                qide = qide.at[:, cm].set(W(mtt, qi[:, cm], qide[:, cm]))
                nide = nide.at[:, cm].set(W(mtt, ni[:, cm], nide[:, cm]))
                diag["dlfm"] = diag["dlfm"].at[:, cm2].set(
                    W(mtt, -du[:, cm2] * qcde[:, cm],
                      diag["dlfm"][:, cm2]))
                diag["dlfn"] = diag["dlfn"].at[:, cm2].set(
                    W(mtt, -du[:, cm2] * ncde[:, cm],
                      diag["dlfn"][:, cm2]))
                diag["difm"] = diag["difm"].at[:, cm2].set(
                    W(mtt, -du[:, cm2] * qide[:, cm],
                      diag["difm"][:, cm2]))
                diag["difn"] = diag["difn"].at[:, cm2].set(
                    W(mtt, -du[:, cm2] * nide[:, cm],
                      diag["difn"][:, cm2]))

            # rain size distribution at k-1 (3033-3059)
            (lamr, n0r, nr, umr, unr) = _rain_sd(
                act, cm, qr, nr, lamr, n0r, umr, unr, arn, rho, mp,
                pre=False)

            # snow at k-1 (3063-3110) + trsp/no-detrain fixups
            (lams, n0s, ns, ums, uns) = _snow_sd(
                act, cm, c, qni, ns, lams, n0s, ums, uns, asn, rho,
                mp, dcs, sp054=True)
            nsde = nsde.at[:, c].set(W(act, ns[:, cm], nsde[:, c]))
            diag["dsfn"] = diag["dsfn"].at[:, cm].set(
                W(act, -du[:, cm] * nsde[:, c], diag["dsfn"][:, cm]))
            # (nsadj is computed but unused in the Fortran)
            if bool(jnp.any(mtt)) and c - 2 >= 0:
                cm2 = c - 2
                qnide = qnide.at[:, cm].set(W(mtt, 0.0, qnide[:, cm]))
                nsde = nsde.at[:, cm].set(W(mtt, 0.0, nsde[:, cm]))
                diag["dsfm"] = diag["dsfm"].at[:, cm2].set(
                    W(mtt, 0.0, diag["dsfm"][:, cm2]))
                diag["dsfn"] = diag["dsfn"].at[:, cm2].set(
                    W(mtt, 0.0, diag["dsfn"][:, cm2]))

            # graupel at k-1 (3113-3143)
            (lamg, n0g, ng, umg, ung) = _graupel_sd(
                act, cm, c, qg, ng, lamg, n0g, umg, ung, agn, rho, mp,
                sp054=True)

            # ---- cleanups below the big if (3149-3181) ----
            (qni, ns, qr, nr, qg, ng, qi, ni, qc, nc) = _cleanup(
                ltrue, cm, qni, ns, qr, nr, qg, ng, qi, ni, qc, nc)

    out = dict(qc=qc, qi=qi, nc=nc, ni=ni, qcde=qcde, qide=qide,
               qnide=qnide, ncde=ncde, nide=nide, nsde=nsde, qni=qni,
               qr=qr, ns=ns, nr=nr, qg=qg, ng=ng, rprd=rprd, sprd=sprd,
               frz=frz, wu=wu, lamc=lamc, pgam=pgam)
    out.update(diag)
    return out


_DIAG_NAMES = ["autolm", "accrlm", "bergnm", "fhtimm", "fhtctm",
               "fhmlm", "hmpim", "accslm", "dlfm", "autoln", "accrln",
               "bergnn", "fhtimn", "fhtctn", "fhmln", "accsln",
               "activn", "dlfn", "autoim", "accsim", "difm", "nuclin",
               "autoin", "accsin", "hmpin", "difn", "trspcm", "trspcn",
               "trspim", "trspin", "accgrm", "accglm", "accgslm",
               "accgsrm", "accgirm", "accgrim", "accgrsm", "accgsln",
               "accgsrn", "accgirn", "accsrim", "acciglm", "accigrm",
               "accsirm", "accigln", "accigrn", "accsirn", "accgln",
               "accgrn", "accilm", "acciln", "fallrm", "fallsm",
               "fallgm", "fallrn", "fallsn", "fallgn", "fhmrm", "dsfm",
               "dsfn"]


def _cleanup(ltrue, cm, qni, ns, qr, nr, qg, ng, qi, ni, qc, nc):
    """Fortran 3149-3181: zero tiny mixing ratios + number floors at
    level k-1 (runs for every k regardless of the plume mask)."""
    m = ltrue & (qni[:, cm] < QSMALL)
    qni = qni.at[:, cm].set(W(m, 0.0, qni[:, cm]))
    ns = ns.at[:, cm].set(W(m, 0.0, ns[:, cm]))
    m = ltrue & (qr[:, cm] < QSMALL)
    qr = qr.at[:, cm].set(W(m, 0.0, qr[:, cm]))
    nr = nr.at[:, cm].set(W(m, 0.0, nr[:, cm]))
    m = ltrue & (qg[:, cm] < QSMALL)
    qg = qg.at[:, cm].set(W(m, 0.0, qg[:, cm]))
    ng = ng.at[:, cm].set(W(m, 0.0, ng[:, cm]))
    m = ltrue & (qi[:, cm] < QSMALL)
    qi = qi.at[:, cm].set(W(m, 0.0, qi[:, cm]))
    ni = ni.at[:, cm].set(W(m, 0.0, ni[:, cm]))
    m = ltrue & (qc[:, cm] < QSMALL)
    qc = qc.at[:, cm].set(W(m, 0.0, qc[:, cm]))
    nc = nc.at[:, cm].set(W(m, 0.0, nc[:, cm]))
    nr = nr.at[:, cm].set(W(ltrue, jnp.maximum(nr[:, cm], 0.0),
                            nr[:, cm]))
    ns = ns.at[:, cm].set(W(ltrue, jnp.maximum(ns[:, cm], 0.0),
                            ns[:, cm]))
    ng = ng.at[:, cm].set(W(ltrue, jnp.maximum(ng[:, cm], 0.0),
                            ng[:, cm]))
    ni = ni.at[:, cm].set(W(ltrue, jnp.maximum(ni[:, cm], 0.0),
                            ni[:, cm]))
    nc = nc.at[:, cm].set(W(ltrue, jnp.maximum(nc[:, cm], 0.0),
                            nc[:, cm]))
    return qni, ns, qr, nr, qg, ng, qi, ni, qc, nc


def _rain_sd(act, idx, qarr, narr, lamr, n0r, umr, unr, arn, rho, mp,
             pre):
    """Rain size distribution + fallspeeds at level `idx` (both the
    pre-process block 1400-1428 on qric/nric and the post-integration
    block 3033-3059 on qr/nr share this exact code)."""
    m = act & (qarr[:, idx] >= QSMALL)
    q_s = W(m, qarr[:, idx], 1.0)
    lam = (PI * RHOW * W(m, narr[:, idx], 0.0) / q_s) ** (1.0 / 3.0)
    lammax = 1.0 / 150.0e-6
    lammin = 1.0 / 3000.0e-6
    lo = m & (lam < lammin)
    hi = m & (lam > lammax)
    lam = W(lo, lammin, W(hi, lammax, lam))
    if pre:
        n0 = W(lo | hi, _divs(lam ** 4.0 * q_s, PI * RHOW),
               W(m, narr[:, idx], 0.0) * lam)
    else:
        n0 = W(lo | hi, _divs(lam ** 4 * q_s, PI * RHOW),
               W(m, narr[:, idx], 0.0) * lam)
    narr = narr.at[:, idx].set(W(lo | hi, n0 / lam, narr[:, idx]))
    lam_s = W(m, lam, 1.0)
    un = jnp.minimum(arn[:, idx] * mp["g_1br"] / lam_s ** BR, 10.0)
    um = jnp.minimum(arn[:, idx] * mp["g_4br"] / (6.0 * lam_s ** BR),
                     10.0)
    lamr = lamr.at[:, idx].set(W(m, lam, W(act, 0.0, lamr[:, idx])))
    n0r = n0r.at[:, idx].set(W(m, n0, W(act, 0.0, n0r[:, idx])))
    umr = umr.at[:, idx].set(W(m, um, W(act, 0.0, umr[:, idx])))
    unr = unr.at[:, idx].set(W(m, un, W(act, 0.0, unr[:, idx])))
    return lamr, n0r, narr, umr, unr


def _snow_sd(act, idx, rho_idx, qarr, narr, lams, n0s, ums, uns, asn,
             rho, mp, dcs, sp054):
    """Snow size distribution + fallspeeds at level `idx` (pre block
    1432-1461 uses rho at the same level and the _r8 0.54 exponent
    with lam**4; the post block 3063-3092 uses rho at level k
    (rho_idx), a single-precision 0.54 and lam**(ds+1)."""
    m = act & (qarr[:, idx] >= QSMALL)
    q_s = W(m, qarr[:, idx], 1.0)
    lam = (mp["g_1ds"] * mp["cs"] * W(m, narr[:, idx], 0.0)
           / q_s) ** (1.0 / mp["ds"])
    lammax = 1.0 / dcs
    lammin = 1.0 / 5000.0e-6
    lo = m & (lam < lammin)
    hi = m & (lam > lammax)
    lam = W(lo, lammin, W(hi, lammax, lam))
    if sp054:
        n0 = W(lo | hi, _divs(lam ** (mp["ds"] + 1.0) * q_s,
                              mp["cs"] * mp["g_1ds"]),
               W(m, narr[:, idx], 0.0) * lam)
        dum = (mp["rhosu"] / rho[:, rho_idx]) ** _EXP_054_SP
    else:
        n0 = W(lo | hi, _divs(lam ** 4.0 * q_s,
                              mp["cs"] * mp["g_1ds"]),
               W(m, narr[:, idx], 0.0) * lam)
        dum = (mp["rhosu"] / rho[:, rho_idx]) ** 0.54
    narr = narr.at[:, idx].set(W(lo | hi, n0 / lam, narr[:, idx]))
    lam_s = W(m, lam, 1.0)
    um = jnp.minimum(asn[:, idx] * mp["g_4bs"] / (6.0 * lam_s ** BS),
                     1.2 * dum)
    un = jnp.minimum(asn[:, idx] * mp["g_1bs"] / lam_s ** BS,
                     1.2 * dum)
    lams = lams.at[:, idx].set(W(m, lam, W(act, 0.0, lams[:, idx])))
    n0s = n0s.at[:, idx].set(W(m, n0, W(act, 0.0, n0s[:, idx])))
    ums = ums.at[:, idx].set(W(m, um, W(act, 0.0, ums[:, idx])))
    uns = uns.at[:, idx].set(W(m, un, W(act, 0.0, uns[:, idx])))
    return lams, n0s, narr, ums, uns


def _graupel_sd(act, idx, rho_idx, qarr, narr, lamg, n0g, umg, ung,
                agn, rho, mp, sp054):
    """Graupel size distribution + fallspeeds at level `idx` (pre
    block 1465-1492 / post block 3113-3143)."""
    m = act & (qarr[:, idx] >= QSMALL)
    q_s = W(m, qarr[:, idx], 1.0)
    lam = (mp["g_1dg"] * mp["cg"] * W(m, narr[:, idx], 0.0)
           / q_s) ** (1.0 / mp["dg"])
    lammax = 1.0 / 20.0e-6
    lammin = 1.0 / 5000.0e-6
    lo = m & (lam < lammin)
    hi = m & (lam > lammax)
    lam = W(lo, lammin, W(hi, lammax, lam))
    n0 = W(lo | hi, _divs(lam ** 4.0 * q_s, mp["g_1dg"] * mp["cg"]),
           W(m, narr[:, idx], 0.0) * lam)
    narr = narr.at[:, idx].set(W(lo | hi, n0 / lam, narr[:, idx]))
    lam_s = W(m, lam, 1.0)
    if sp054:
        dum = (mp["rhosu"] / rho[:, rho_idx]) ** _EXP_054_SP
    else:
        dum = (mp["rhosu"] / rho[:, rho_idx]) ** 0.54
    um = jnp.minimum(agn[:, idx] * mp["g_4bg"] / (6.0 * lam_s ** BG),
                     20.0 * dum)
    un = jnp.minimum(agn[:, idx] * mp["g_1bg"] / lam_s ** BG,
                     20.0 * dum)
    lamg = lamg.at[:, idx].set(W(m, lam, W(act, 0.0, lamg[:, idx])))
    n0g = n0g.at[:, idx].set(W(m, n0, W(act, 0.0, n0g[:, idx])))
    umg = umg.at[:, idx].set(W(m, um, W(act, 0.0, umg[:, idx])))
    ung = ung.at[:, idx].set(W(m, un, W(act, 0.0, ung[:, idx])))
    return lamg, n0g, narr, umg, ung


# ---------------------------------------------------------------------------
# zm_microphysics_adjust
# ---------------------------------------------------------------------------
def zm_microphysics_adjust(jt, msg, delt, dp, qv, dl, dsdt, dqdt,
                           rprd, sprd, dnlf, dif, dnif, dsf, dnsf, zc):
    """zm_microphysics_adjust: evaporate precipitation / detrained
    condensate to prevent negative water vapor. jt 0-based; msg keeps
    its Fortran VALUE. Returns updated (dl, dsdt, dqdt, rprd, sprd,
    dnlf, dif, dnif, dsf, dnsf)."""
    dl = jnp.asarray(dl, dtype=jnp.float64)
    dsdt = jnp.asarray(dsdt, dtype=jnp.float64)
    dqdt = jnp.asarray(dqdt, dtype=jnp.float64)
    rprd = jnp.asarray(rprd, dtype=jnp.float64)
    sprd = jnp.asarray(sprd, dtype=jnp.float64)
    n, pver = dl.shape
    lv_cp = zc["latvap"] / zc["cpair"]
    li_cp = zc["latice"] / zc["cpair"]
    lvi_cp = (zc["latvap"] + zc["latice"]) / zc["cpair"]
    jtmin = int(jnp.min(jt))

    for kf in range(msg + 1, pver + 1):     # Fortran k
        c = kf - 1
        start = (dqdt[:, c] * delt + qv[:, c]) < 0.0
        negadq = W(start, dqdt[:, c] + _divs(qv[:, c], delt), 0.0)
        dqdt = dqdt.at[:, c].set(W(start, dqdt[:, c] - negadq,
                                   dqdt[:, c]))
        dpk = dp[:, c]
        for lf in range(kf, jtmin, -1):     # Fortran kk
            kkc = lf - 1
            m0 = start & (kkc >= jt) & (negadq < 0.0)
            dpkk = dp[:, kkc]
            enough = m0 & (rprd[:, kkc] > -negadq * dpk / dpkk)
            # --- precipitation branch ---
            dsdt = dsdt.at[:, c].set(
                W(enough, dsdt[:, c] + negadq * lv_cp, dsdt[:, c]))
            has_rain = enough & (rprd[:, kkc] > sprd[:, kkc])
            rain_short = has_rain & ((rprd[:, kkc] - sprd[:, kkc])
                                     < -negadq * dpk / dpkk)
            dsdt = dsdt.at[:, c].set(
                W(rain_short, dsdt[:, c]
                  + (negadq + (rprd[:, kkc] - sprd[:, kkc])
                     * dpkk / dpk) * li_cp, dsdt[:, c]))
            sprd = sprd.at[:, kkc].set(
                W(rain_short, negadq * dpk / dpkk + rprd[:, kkc],
                  sprd[:, kkc]))
            no_rain = enough & ~has_rain
            sprd = sprd.at[:, kkc].set(
                W(no_rain, sprd[:, kkc] + negadq * dpk / dpkk,
                  sprd[:, kkc]))
            dsdt = dsdt.at[:, c].set(
                W(no_rain, dsdt[:, c] + negadq * li_cp, dsdt[:, c]))
            rprd = rprd.at[:, kkc].set(
                W(enough, rprd[:, kkc] + negadq * dpk / dpkk,
                  rprd[:, kkc]))
            not_enough = m0 & ~enough
            negadq = W(enough, 0.0,
                       W(not_enough,
                         rprd[:, kkc] * dpkk / dpk + negadq, negadq))
            dsdt = dsdt.at[:, c].set(
                W(not_enough, dsdt[:, c]
                  - rprd[:, kkc] * lv_cp * dpkk / dpk
                  - sprd[:, kkc] * li_cp * dpkk / dpk, dsdt[:, c]))
            sprd = sprd.at[:, kkc].set(W(not_enough, 0.0, sprd[:, kkc]))
            rprd = rprd.at[:, kkc].set(W(not_enough, 0.0, rprd[:, kkc]))
            # --- detrained condensate branch ---
            m1 = m0 & (negadq < 0.0)
            dl_enough = m1 & (dl[:, kkc] > -negadq * dpk / dpkk)
            dsdt = dsdt.at[:, c].set(
                W(dl_enough, dsdt[:, c] + negadq * lv_cp, dsdt[:, c]))
            dl_s = W(dl_enough, dl[:, kkc], 1.0)
            dnlf = dnlf.at[:, kkc].set(
                W(dl_enough, dnlf[:, kkc]
                  * (1.0 + negadq * dpk / dpkk / dl_s), dnlf[:, kkc]))
            dl = dl.at[:, kkc].set(
                W(dl_enough, dl[:, kkc] + negadq * dpk / dpkk,
                  dl[:, kkc]))
            m2 = m1 & ~dl_enough
            negadq = W(dl_enough, 0.0,
                       W(m2, negadq + dl[:, kkc] * dpkk / dpk, negadq))
            dsdt = dsdt.at[:, c].set(
                W(m2, dsdt[:, c] - dl[:, kkc] * dpkk / dpk * lv_cp,
                  dsdt[:, c]))
            dl = dl.at[:, kkc].set(W(m2, 0.0, dl[:, kkc]))
            dnlf = dnlf.at[:, kkc].set(W(m2, 0.0, dnlf[:, kkc]))
            dif_enough = m2 & (dif[:, kkc] > -negadq * dpk / dpkk)
            dsdt = dsdt.at[:, c].set(
                W(dif_enough, dsdt[:, c] + negadq * lvi_cp,
                  dsdt[:, c]))
            dif_s = W(dif_enough, dif[:, kkc], 1.0)
            dnif = dnif.at[:, kkc].set(
                W(dif_enough, dnif[:, kkc]
                  * (1.0 + negadq * dpk / dpkk / dif_s), dnif[:, kkc]))
            dif = dif.at[:, kkc].set(
                W(dif_enough, dif[:, kkc] + negadq * dpk / dpkk,
                  dif[:, kkc]))
            m3 = m2 & ~dif_enough
            negadq = W(dif_enough, 0.0,
                       W(m3, negadq + dif[:, kkc] * dpkk / dpk,
                         negadq))
            dsdt = dsdt.at[:, c].set(
                W(m3, dsdt[:, c] - dif[:, kkc] * dpkk / dpk * lvi_cp,
                  dsdt[:, c]))
            dif = dif.at[:, kkc].set(W(m3, 0.0, dif[:, kkc]))
            dnif = dnif.at[:, kkc].set(W(m3, 0.0, dnif[:, kkc]))
            dsf_enough = m3 & (dsf[:, kkc] > -negadq * dpk / dpkk)
            dsdt = dsdt.at[:, c].set(
                W(dsf_enough, dsdt[:, c] + negadq * lvi_cp,
                  dsdt[:, c]))
            dsf_s = W(dsf_enough, dsf[:, kkc], 1.0)
            dnsf = dnsf.at[:, kkc].set(
                W(dsf_enough, dnsf[:, kkc]
                  * (1.0 + negadq * dpk / dpkk / dsf_s), dnsf[:, kkc]))
            dsf = dsf.at[:, kkc].set(
                W(dsf_enough, dsf[:, kkc] + negadq * dpk / dpkk,
                  dsf[:, kkc]))
            m4 = m3 & ~dsf_enough
            negadq = W(dsf_enough, 0.0,
                       W(m4, negadq + dsf[:, kkc] * dpkk / dpk,
                         negadq))
            dsdt = dsdt.at[:, c].set(
                W(m4, dsdt[:, c] - dsf[:, kkc] * dpkk / dpk * lvi_cp,
                  dsdt[:, c]))
            dsf = dsf.at[:, kkc].set(W(m4, 0.0, dsf[:, kkc]))
            dnsf = dnsf.at[:, kkc].set(W(m4, 0.0, dnsf[:, kkc]))
        dqdt = dqdt.at[:, c].set(
            W(start & (negadq < 0.0), dqdt[:, c] - negadq, dqdt[:, c]))
    return dl, dsdt, dqdt, rprd, sprd, dnlf, dif, dnif, dsf, dnsf
