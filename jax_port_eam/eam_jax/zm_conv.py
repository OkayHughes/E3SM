"""Port of the ZM deep-convection main routine,
eam/src/physics/cam/zm/zm_conv.F90 zm_conv_main: dilute CAPE + DCAPE
trigger, active-column gather, updraft/downdraft plume properties
(zm_cloud_properties = zm_calc_fractional_entrainment + plume loops +
zm_downdraft_properties), CAPE closure (zm_closure), cloud-base
mass-flux limiting/scaling, output tendencies (zm_calc_output_tend),
and the precip / reserved-liquid integrals.

PORT_NOTES
----------
- SCOPE: zm_param%zm_microp = False only (the simple in-plume
  condensation path; itnum = 1). The EAMv3 default is
  zmconv_microp=.true., but zm_microphysics.F90 (3300 lines) needs the
  aerosol activation stack and is deferred (PORTING_PLAN.md row 9).
  The Fortran golden harness runs with the same flag false, so goldens
  and port match. MCSP and zm_aero paths are inactive by construction
  (mcsp lives in zm_conv_intr/zm_conv_mcsp; aero is only touched under
  zm_microp).
- zm_conv_evap (below-cloud Sundqvist evaporation of convective
  precip + snow melt/production) IS ported here, together with the
  real cloud_fraction.F90 cldfrc_fice it calls (T-based ice/snow
  partition; top_lev effectively 1 -- in production top_lev =
  trop_cloud_top_lev sits above the 40 hPa limcnv cap where all ZM
  precip fluxes are identically zero, so the choice cannot affect
  outputs). Both zm_param%old_snow branches are ported and goldened;
  with zm_microp=False the not-old_snow branch has prdsnow = 0
  (microp_st%sprd needs zm_microphysics, PORTING_PLAN.md row 9), so
  its snow flux is identically zero and its final flxsnow<=flxprec
  protection loop is unreachable (ported verbatim, validated only for
  no-op). tend_s/tend_q are declared inout in the Fortran but every
  (i,k) is overwritten before being read, so the port treats them as
  pure outputs. zm_param%ke is the single zmconv_ke evaporation
  efficiency (no zmconv_ke_lnd exists in EAMv3 zm_conv.F90).
  pergro_active = False (no PERGRO).
- Index convention: 0-based, level 0 = model top. Fortran midpoint k
  maps to python k-1. msg = limcnv-1 keeps its Fortran VALUE (it is
  the count of excluded top levels; Fortran loops k = msg+1..pver are
  python j = msg..pver-1). The Fortran 1-based initial plume top
  jt = max(lel, limcnv+1) becomes jt0 = max(lel0, msg+1).
- Gathering: active columns (cape/dcape trigger) are selected with
  concrete numpy indices; plume/closure/tendency math runs on the
  gathered rows only, then results are scattered exactly as the
  Fortran does. GATHERED outputs (msemax_klev_g, jt, mflx_up, entr_up,
  detr_up, mflx_dn, entr_dn, p_del, dsubcld, ql) are returned in
  gathered row order, with the same fill values the Fortran leaves in
  rows >= lengath (zeros; jt rows = pver-1 from its whole-array init;
  msemax_klev_g rows = 0 i.e. Fortran 1). SCATTERED outputs (jctop,
  jcbot, prec, heat, qtnd, cape, dcape, mcon, pflx, zdu, rprd, dlf,
  rliq) are in column order.
- The CAPE step reuses eam_jax.zm_cape.compute_dilute_cape verbatim
  (call 1: calc_msemax_klev=True; DCAPE call 2 on t_star/q_star with
  the frozen launch level), including its zm_const (zvir = 1.608
  hardcoded in zm_conv_types).
- Sequential level sweeps (updraft MSE, LCL search, in-plume
  condensate ql, downdraft h/s/q recursions, dsdt/dqdt below cloud
  base, prec integral) are python loops in the exact Fortran order
  with per-column masks; where a level's update reads the level below
  (k+1), the loop direction preserves the data flow. Independent-in-k
  loops are vectorized. All reductions accumulate in Fortran loop
  order (no jnp.sum where order matters).
- Faithful oddities preserved:
  * the "nonsensical" -alfa*lambda_max/lambda_max downdraft seed
    (identical to -alfa in IEEE for lambda_max > 0; kept simple);
  * h_dnd at the downdraft source level jd is seeded with
    h_env(jd-1) and then OVERWRITTEN by the k-ascending recursion,
    which reproduces h_env(jd-1) only up to roundoff
    (dz*(m/dz) != m) - the recursion result is what's kept;
  * downdraft recursions divide by min(mflx_dn, -1e-20), producing
    O(1e24) garbage in h_dnd/s_dnd outside [jd, jb]; those entries are
    never read (float64 carries them without overflow) but s_dnd IS
    fed to zm_closure/zm_calc_output_tend multiplied by mflx_dn == 0
    exactly, so 0 * garbage = 0 keeps results finite and exact;
  * j0 (detrainment start) is intent(out) and only defined when the
    MSE-min search fires; for gathered columns jt < jb always holds
    (cape > 0 puts the equilibrium level above the launch level), so
    the first visited level always fires - the port asserts this
    instead of reproducing undefined behaviour;
  * zm_closure's dry/moist buoyancy-derivative formulas use literal
    0.608/1.608 (not zvir), and theta factors 1000/p_mid in hPa;
  * beta = 0 in zm_closure (compile-time parameter) - terms kept as
    beta*ql(k) + (1-beta)*ql(k+1) with beta = 0.0 verbatim.
- Units follow the Fortran: pressures converted to hPa (mb) on entry,
  p_del_in stays Pa inside the precip/rliq integrals, mass fluxes are
  normalized by cloud-base mass flux until the closure scaling, and
  entr/detr/cu/rprd/evp are converted from 1/m to 1/mb before the
  closure. heat = dsdt*cpair [W/kg]; prec/rliq in m/s.
- Everything float64, eager-mode (concrete gather indices and python
  min/max over per-chunk index arrays; not jit-compatible as written).
"""

import numpy as np
import jax.numpy as jnp

from .wv_sat import qsat as _qsat_pa
from .zm_cape import compute_dilute_cape, make_zm_const, qsat_hpa
from .zm_cape import make_zm_param as _make_cape_param

__all__ = ["make_zm_const", "make_zm_param", "zm_conv_main",
           "cldfrc_fice", "zm_conv_evap"]

# zm_conv module parameters (zm_conv.F90)
CAPE_THRESHOLD_OLD = 70.0     # [J/kg] pre-DCAPE trigger threshold
CAPE_THRESHOLD_NEW = 0.0      # [J/kg] DCAPE-mode cape threshold
DCAPE_THRESHOLD = 0.0         # [J/kg/s] dcape trigger threshold
INTERP_DIFF_MIN = 1.0e-6      # interface log-interpolation threshold
OMSM = 0.99999                # round-off protection factor
SMALL = 1.0e-20               # mass-flux normalization floor
# zm_cloud_properties parameters
MU_MIN = 0.02                 # minimum normalized updraft mass flux
HU_DIFF_MIN = -2000.0         # limit on hu-hsthat at cloud top
# zm_calc_fractional_entrainment parameters
LAMBDA_LIMIT_MIN = 0.0
LAMBDA_LIMIT_MAX = 0.0002
LAMBDA_THRESHOLD = 1.0e-6
# zm_closure parameter
BETA = 0.0                    # proportion of liquid from layer below


def make_zm_param(tau=3600.0, alfa=0.14, ke=2.5e-6, dmpdz=-0.7e-3,
                  tpert_fix=True, tpert_fac=2.0, tiedke_add=0.8,
                  c0_lnd=0.0020, c0_ocn=0.0020, num_cin=1, limcnv=1,
                  mx_bot_lyr_adj=1, trig_dcape=True, trig_ull=True,
                  clos_dyn_adj=True, no_deep_pbl=False,
                  old_snow=False):
    """zm_param_t subset used by zm_conv_main with zm_microp=False.
    Defaults are the EAMv3 phys="default" namelist values (see the
    golden generator header); limcnv is the 1-based interface index
    from the 40 hPa rule and must be supplied per grid."""
    zp = _make_cape_param(dmpdz=dmpdz, tiedke_add=tiedke_add,
                          tpert_fac=tpert_fac,
                          mx_bot_lyr_adj=mx_bot_lyr_adj,
                          num_cin=num_cin, tpert_fix=tpert_fix,
                          trig_ull=trig_ull, trig_dcape=trig_dcape)
    zp.update(tau=float(tau), alfa=float(alfa), ke=float(ke),
              c0_lnd=float(c0_lnd), c0_ocn=float(c0_ocn),
              limcnv=int(limcnv), clos_dyn_adj=bool(clos_dyn_adj),
              no_deep_pbl=bool(no_deep_pbl), old_snow=bool(old_snow),
              zm_microp=False)
    return zp


def _interface_interp(mid, msg):
    """Interfacial (q,s) values (zm_conv_main): log interpolation when
    the relative jump exceeds INTERP_DIFF_MIN, else arithmetic mean.
    Interface j (0-based) sits between midpoints j-1 and j; entries
    j <= msg keep the midpoint value (Fortran loop is k = msg+2..pver
    and initialization is the midpoint copy)."""
    n, pver = mid.shape
    out = jnp.array(mid)
    a = mid[:, :-1]   # level k-1
    b = mid[:, 1:]    # level k
    difr = jnp.where((a > 0.0) | (b > 0.0),
                     jnp.abs((b - a) / jnp.maximum(a, b)), 0.0)
    safe = jnp.abs(a - b) > 0.0
    logint = jnp.where(safe,
                       jnp.log(jnp.where(safe, a / b, 1.0)) * a * b
                       / jnp.where(safe, a - b, 1.0), b)
    vals = jnp.where(difr > INTERP_DIFF_MIN, logint, 0.5 * (a + b))
    ki = jnp.arange(1, pver)[None, :]
    out = out.at[:, 1:].set(jnp.where(ki >= msg + 1, vals, out[:, 1:]))
    return out


def zm_calc_fractional_entrainment(jb, jt0, j0, z_mid, z_int, dz,
                                   h_env, h_env_sat, msg, zc):
    """Fractional entrainment lambda(z) via the 4-term Taylor series of
    ZM95 eq (A6), with the detrainment-level adjustment. Returns
    (j0, lambda, lambda_max)."""
    n, pver = h_env.shape
    rows = jnp.arange(n)
    k1 = jnp.zeros((n, pver))
    i2 = jnp.zeros((n, pver))
    i3 = jnp.zeros((n, pver))
    i4 = jnp.zeros((n, pver))
    lam_tmp = jnp.zeros((n, pver))
    lam = jnp.zeros((n, pver))
    h_env_jb = h_env[rows, jb]

    # Taylor-series integrals (sequential, bottom -> top)
    for j in range(pver - 2, msg - 1, -1):
        m = (j < jb) & (j >= jt0)
        k1j = k1[:, j + 1] + (h_env_jb - h_env[:, j]) * dz[:, j]
        k1 = k1.at[:, j].set(jnp.where(m, k1j, k1[:, j]))
        ihat = 0.5 * (k1[:, j + 1] + k1[:, j])
        i2j = i2[:, j + 1] + ihat * dz[:, j]
        i2 = i2.at[:, j].set(jnp.where(m, i2j, i2[:, j]))
        idag = 0.5 * (i2[:, j + 1] + i2[:, j])
        i3j = i3[:, j + 1] + idag * dz[:, j]
        i3 = i3.at[:, j].set(jnp.where(m, i3j, i3[:, j]))
        iprm = 0.5 * (i3[:, j + 1] + i3[:, j])
        i4j = i4[:, j + 1] + iprm * dz[:, j]
        i4 = i4.at[:, j].set(jnp.where(m, i4j, i4[:, j]))

    # re-initialized minimum MSE over [j0, jb]
    h_env_min = jnp.full(n, 1.0e6)
    for j in range(msg, pver):
        m = (j >= j0) & (j <= jb) & (h_env[:, j] <= h_env_min)
        h_env_min = jnp.where(m, h_env[:, j], h_env_min)

    # approximate lambda(z); k1 zeroing outside the plume matters for
    # the branch condition k1 > expnum*dz
    for j in range(msg + 1, pver):
        outside = (j < jt0) | (j >= jb)
        k1 = k1.at[:, j].set(jnp.where(outside, 0.0, k1[:, j]))
        expnum = jnp.where(
            outside, 0.0,
            h_env_jb - (h_env_sat[:, j - 1]
                        * (z_int[:, j] - z_mid[:, j])
                        + h_env_sat[:, j]
                        * (z_mid[:, j - 1] - z_int[:, j]))
            / (z_mid[:, j - 1] - z_mid[:, j]))
        m = (((h_env_jb - h_env_min) > 100.0) & (expnum > 0.0)
             & (k1[:, j] > expnum * dz[:, j]))
        k1j = jnp.where(m, k1[:, j], 1.0)      # safe denominator
        tmp = expnum / k1j
        lam_j = (tmp
                 + i2[:, j] / k1j * tmp**2
                 + (2.0 * i2[:, j]**2 - k1j * i3[:, j]) / k1j**2
                 * tmp**3
                 + (-5.0 * k1j * i2[:, j] * i3[:, j]
                    + 5.0 * i2[:, j]**3 + k1j**2 * i4[:, j])
                 / k1j**3 * tmp**4)
        lam_j = jnp.clip(lam_j, LAMBDA_LIMIT_MIN, LAMBDA_LIMIT_MAX)
        lam_tmp = lam_tmp.at[:, j].set(
            jnp.where(m, lam_j, lam_tmp[:, j]))

    # move detrainment level down if entrainment is too weak
    lt_j0 = lam_tmp[rows, j0]
    lt_j0p1 = lam_tmp[rows, jnp.minimum(j0 + 1, pver - 1)]
    move = ((j0 < jb) & (lt_j0 < LAMBDA_THRESHOLD) & (lt_j0p1 > lt_j0))
    j0 = jnp.where(move, j0 + 1, j0)

    # entrainment must not increase above the detrainment level
    for j in range(msg + 1, pver):
        m = (j >= jt0) & (j <= j0)
        lam_tmp = lam_tmp.at[:, j].set(
            jnp.where(m, jnp.maximum(lam_tmp[:, j], lam_tmp[:, j - 1]),
                      lam_tmp[:, j]))

    lambda_max = lam_tmp[rows, j0]
    lam = lam.at[rows, jb].set(lambda_max)
    ki = jnp.arange(pver)[None, :]
    lam = jnp.where((ki >= j0[:, None]) & (ki <= jb[:, None]),
                    lam_tmp[rows, j0][:, None], lam)
    lam = jnp.where((ki < j0[:, None]) & (ki >= jt0[:, None]),
                    lam_tmp, lam)
    return j0, lam, lambda_max


def zm_downdraft_properties(jb, jt0, j0, z_int, dz, s_mid, q_mid,
                            h_env, lambda_max, qsthat, hsthat, gamhat,
                            rprd, mflx_up, msg, zc, zp,
                            s_dnd, q_dnd, q_dnd_sat, h_dnd):
    """Downdraft mass flux, entrainment, MSE/q/s recursions and
    evaporation. Returns (jt0, jd, mflx_dn, entr_dn, s_dnd, q_dnd,
    h_dnd, q_dnd_sat, evp, totevp)."""
    n, pver = s_mid.shape
    rows = jnp.arange(n)
    lm_pos = lambda_max > 0.0
    mflx_dn = jnp.zeros((n, pver))
    entr_dn = jnp.zeros((n, pver))
    evp = jnp.zeros((n, pver))
    totevp = jnp.zeros(n)

    jt0 = jnp.minimum(jt0, jb - 1)
    jd = jnp.maximum(j0, jt0 + 1)
    jd = jnp.minimum(jd, jb)
    h_dnd = h_dnd.at[rows, jd].set(h_env[rows, jd - 1])
    seed = (jd < jb) & lm_pos
    lm_safe = jnp.where(lm_pos, lambda_max, 1.0)
    mflx_dn = mflx_dn.at[rows, jd].set(
        jnp.where(seed, -zp["alfa"] * lm_safe / lm_safe,
                  mflx_dn[rows, jd]))

    # mass-flux profile below the source level
    ki = jnp.arange(pver)[None, :]
    m = (ki > jd[:, None]) & (ki <= jb[:, None]) & lm_pos[:, None]
    dz_tmp = z_int[rows, jd][:, None] - z_int[:, :pver]
    dz_safe = jnp.where(m, dz_tmp, 1.0)
    prof = (-zp["alfa"] / (2.0 * lm_safe[:, None])
            * (jnp.exp(2.0 * lm_safe[:, None] * dz_safe) - 1.0)
            / dz_safe)
    mflx_dn = jnp.where(m, prof, mflx_dn)

    # scale so the net cloud-base flux is not negative; ratmjb is
    # computed from the unscaled base value for every level incl. jb
    scale_m = lm_pos & (jd < jb)
    denom = jnp.where(scale_m, mflx_dn[rows, jb], -1.0)
    ratmjb = jnp.minimum(jnp.abs(mflx_up[rows, jb] / denom), 1.0)
    m = ((ki >= jt0[:, None]) & (ki <= jb[:, None])
         & scale_m[:, None])
    mflx_dn = jnp.where(m, mflx_dn * ratmjb[:, None], mflx_dn)

    # entrainment + MSE recursion (sequential, top -> bottom); the
    # min(mflx_dn, -1e-20) floor makes h_dnd garbage outside [jd, jb]
    # exactly as in the Fortran (never read there)
    for j in range(msg, pver):
        m = (j >= jt0) & lm_pos
        ed_new = (mflx_dn[:, j - 1] - mflx_dn[:, j]) / dz[:, j - 1]
        entr_dn = entr_dn.at[:, j - 1].set(
            jnp.where(m, ed_new, entr_dn[:, j - 1]))
        mdt = jnp.minimum(mflx_dn[:, j], -SMALL)
        hd_new = (mflx_dn[:, j - 1] * h_dnd[:, j - 1]
                  - dz[:, j - 1] * entr_dn[:, j - 1] * h_env[:, j - 1]
                  ) / mdt
        h_dnd = h_dnd.at[:, j].set(jnp.where(m, hd_new, h_dnd[:, j]))

    # downdraft saturation specific humidity
    m = ((ki >= jd[:, None]) & (ki <= jb[:, None]) & lm_pos[:, None]
         & (jd < jb)[:, None] & (ki >= msg + 1))
    qds = qsthat + gamhat * (h_dnd - hsthat) \
        / (zc["latvap"] * (1.0 + gamhat))
    q_dnd_sat = jnp.where(m, qds, q_dnd_sat)

    # source-level values (unconditional, as in the Fortran)
    q_dnd = q_dnd.at[rows, jd].set(q_dnd_sat[rows, jd])
    s_dnd = s_dnd.at[rows, jd].set(
        (h_dnd[rows, jd] - zc["latvap"] * q_dnd[rows, jd])
        / zc["cpair"])

    # evaporation + s recursion (sequential, top -> bottom)
    for j in range(msg + 1, pver - 1):
        m = (j >= jd) & (j < jb) & lm_pos
        q_dnd = q_dnd.at[:, j + 1].set(
            jnp.where(m, q_dnd_sat[:, j + 1], q_dnd[:, j + 1]))
        evp_new = (-entr_dn[:, j] * q_mid[:, j]
                   + (mflx_dn[:, j] * q_dnd[:, j]
                      - mflx_dn[:, j + 1] * q_dnd[:, j + 1])
                   / dz[:, j])
        evp_new = jnp.maximum(evp_new, 0.0)
        mdt = jnp.minimum(mflx_dn[:, j + 1], -SMALL)
        sd_new = ((zc["latvap"] / zc["cpair"] * evp_new
                   - entr_dn[:, j] * s_mid[:, j]) * dz[:, j]
                  + mflx_dn[:, j] * s_dnd[:, j]) / mdt
        s_dnd = s_dnd.at[:, j + 1].set(
            jnp.where(m, sd_new, s_dnd[:, j + 1]))
        evp = evp.at[:, j].set(jnp.where(m, evp_new, evp[:, j]))
        totevp = totevp - jnp.where(
            m, dz[:, j] * entr_dn[:, j] * q_mid[:, j], 0.0)

    totevp = (totevp + mflx_dn[rows, jd] * q_dnd[rows, jd]
              - mflx_dn[rows, jb] * q_dnd[rows, jb])
    return (jt0, jd, mflx_dn, entr_dn, s_dnd, q_dnd, h_dnd,
            q_dnd_sat, evp, totevp)


def zm_cloud_properties(p_mid, z_mid, z_int, t_mid, s_mid, s_int,
                        q_mid, landfrac, tpert_g, jb, lel, msg, zc,
                        zp):
    """Updraft/downdraft plume properties (zm_microp=False branch,
    itnum=1). All inputs gathered, 0-based indices. Returns a dict."""
    n, pver = t_mid.shape
    rows = jnp.arange(n)
    ki = jnp.arange(pver)[None, :]

    c0mask = (zp["c0_ocn"] * (1.0 - landfrac)
              + zp["c0_lnd"] * landfrac)
    dz = z_int[:, :pver] - z_int[:, 1:]

    est, qst = qsat_hpa(t_mid, p_mid)
    qst = jnp.where(p_mid - est <= 0.0, 1.0, qst)
    gamma = (qst * (1.0 + qst / zc["epsilo"]) * zc["epsilo"]
             * zc["latvap"] / (zc["rdair"] * t_mid**2)
             * zc["latvap"] / zc["cpair"])

    s_upd = jnp.array(s_mid)
    s_dnd = jnp.array(s_mid)
    q_dnd = jnp.array(q_mid)
    q_dnd_sat = jnp.array(q_mid)
    q_upd = jnp.array(q_mid)
    h_env = (zc["cpair"] * t_mid + zc["grav"] * z_mid
             + zc["latvap"] * q_mid)
    h_env_sat = (zc["cpair"] * t_mid + zc["grav"] * z_mid
                 + zc["latvap"] * qst)
    h_upd = jnp.array(h_env)
    h_dnd = jnp.array(h_env)

    # midpoint -> interface interpolation for qst/hst/gamma
    a, b = qst[:, :-1], qst[:, 1:]
    qsthat_v = jnp.where(jnp.abs(a - b) > INTERP_DIFF_MIN,
                         jnp.log(a / b) * a * b
                         / jnp.where(a == b, 1.0, a - b), b)
    qsthat = jnp.array(qst)
    qsthat = qsthat.at[:, 1:].set(
        jnp.where(jnp.arange(1, pver)[None, :] >= msg + 1,
                  qsthat_v, qsthat[:, 1:]))
    hsthat = jnp.array(h_env_sat)
    hsthat = hsthat.at[:, 1:].set(
        jnp.where(jnp.arange(1, pver)[None, :] >= msg + 1,
                  zc["cpair"] * s_int[:, 1:]
                  + zc["latvap"] * qsthat[:, 1:], hsthat[:, 1:]))
    a, b = gamma[:, :-1], gamma[:, 1:]
    gamhat_v = jnp.where(jnp.abs(a - b) > INTERP_DIFF_MIN,
                         jnp.log(a / b) * a * b
                         / jnp.where(a == b, 1.0, a - b), b)
    gamhat = jnp.array(gamma)
    gamhat = gamhat.at[:, 1:].set(
        jnp.where(jnp.arange(1, pver)[None, :] >= msg + 1,
                  gamhat_v, gamhat[:, 1:]))

    # initial plume top / downdraft level / LCL
    jt0 = jnp.minimum(jnp.maximum(lel, msg + 1), pver - 1)
    jlcl = jnp.array(lel)

    # level of minimum saturated MSE -> detrainment start j0
    h_env_min = jnp.full(n, 1.0e6)
    j0 = jnp.zeros(n, dtype=jnp.int64)
    found = jnp.zeros(n, dtype=bool)
    for j in range(msg, pver):
        m = ((h_env_sat[:, j] <= h_env_min) & (j >= jt0) & (j <= jb))
        h_env_min = jnp.where(m, h_env_sat[:, j], h_env_min)
        j0 = jnp.where(m, j, j0)
        found = found | m
    if not bool(jnp.all(found)):
        raise RuntimeError(
            "zm_cloud_properties: j0 undefined for a gathered column "
            "(jt > jb); the Fortran reads uninitialized memory here")
    j0 = jnp.minimum(j0, jb - 2)
    j0 = jnp.maximum(j0, jt0 + 2)
    j0 = jnp.minimum(j0, pver - 1)

    # updraft MSE/DSE seed over [jt, jb]
    plume = (ki >= jt0[:, None]) & (ki <= jb[:, None])
    if zp["tpert_fix"]:
        dh = zc["cpair"] * zp["tiedke_add"]
        ds = zp["tiedke_add"]
        h_upd = jnp.where(plume, h_env[rows, jb][:, None] + dh, h_upd)
        s_upd = jnp.where(plume, s_mid[rows, jb][:, None] + ds, s_upd)
    else:
        pert = zp["tiedke_add"] + zp["tpert_fac"] * tpert_g
        h_upd = jnp.where(plume, (h_env[rows, jb]
                                  + zc["cpair"] * pert)[:, None],
                          h_upd)
        s_upd = jnp.where(plume, (s_mid[rows, jb] + pert)[:, None],
                          s_upd)

    j0, lam, lambda_max = zm_calc_fractional_entrainment(
        jb, jt0, j0, z_mid, z_int, dz, h_env, h_env_sat, msg, zc)
    lm_pos = lambda_max > 0.0
    lm_safe = jnp.where(lm_pos, lambda_max, 1.0)

    # ---- single-iteration (itnum=1) plume property pass ----
    cu = jnp.zeros((n, pver))
    ql = jnp.zeros((n, pver))
    totpcp = jnp.zeros(n)
    mflx_up = jnp.zeros((n, pver))
    entr_up = jnp.zeros((n, pver))
    detr_up = jnp.zeros((n, pver))

    mflx_up = mflx_up.at[rows, jb].set(jnp.where(lm_pos, 1.0, 0.0))
    entr_up = entr_up.at[rows, jb].set(
        jnp.where(lm_pos, 1.0 / dz[rows, jb], 0.0))
    z_int_jb = z_int[rows, jb]
    for j in range(pver - 1, msg - 1, -1):
        m = lm_pos & (j >= jt0) & (j < jb)
        zuef = z_int[:, j] - z_int_jb
        zuef_s = jnp.where(m, zuef, 1.0)
        rmue = ((1.0 / lm_safe)
                * (jnp.exp(lam[:, j + 1] * zuef_s) - 1.0) / zuef_s)
        mf = ((1.0 / lm_safe)
              * (jnp.exp(lam[:, j] * zuef_s) - 1.0) / zuef_s)
        entr_up = entr_up.at[:, j].set(
            jnp.where(m, (rmue - mflx_up[:, j + 1]) / dz[:, j],
                      entr_up[:, j]))
        detr_up = detr_up.at[:, j].set(
            jnp.where(m, (rmue - mf) / dz[:, j], detr_up[:, j]))
        mflx_up = mflx_up.at[:, j].set(
            jnp.where(m, mf, mflx_up[:, j]))

    kh0 = int(jnp.min(lel))   # khighest (0-based)
    kl0 = int(jnp.max(jb))    # klowest  (0-based)

    # updraft MSE recursion with weak-plume pruning (bottom -> top)
    for j in range(kl0 - 1, kh0 - 1, -1):
        act = (j <= jb - 1) & (j >= lel) & lm_pos
        weak = act & (mflx_up[:, j] < 0.02)
        strong = act & ~weak
        mf_s = jnp.where(mflx_up[:, j] != 0.0, mflx_up[:, j], 1.0)
        hu_new = (mflx_up[:, j + 1] / mf_s * h_upd[:, j + 1]
                  + dz[:, j] / mf_s
                  * (entr_up[:, j] * h_env[:, j]
                     - detr_up[:, j] * h_env_sat[:, j]))
        h_upd = h_upd.at[:, j].set(
            jnp.where(weak, h_env[:, j],
                      jnp.where(strong, hu_new, h_upd[:, j])))
        detr_up = detr_up.at[:, j].set(
            jnp.where(weak, mflx_up[:, j + 1] / dz[:, j],
                      detr_up[:, j]))
        mflx_up = mflx_up.at[:, j].set(
            jnp.where(weak, 0.0, mflx_up[:, j]))
        entr_up = entr_up.at[:, j].set(
            jnp.where(weak, 0.0, entr_up[:, j]))

    # cloud-top search (bottom -> top; tot_frz = 0 with microp off)
    doit = jnp.ones(n, dtype=bool)
    tot_frz = jnp.zeros(n)
    h_upd_jb = h_upd[rows, jb]
    for j in range(kl0 - 2, kh0 - 2, -1):
        if j < 0:
            break
        cond = doit & (j <= jb - 2) & (j >= lel - 1)
        b1 = (cond & (h_upd[:, j] <= hsthat[:, j])
              & (h_upd[:, j + 1] > hsthat[:, j + 1])
              & (mflx_up[:, j] >= MU_MIN))
        b1_low = b1 & (h_upd[:, j] - hsthat[:, j] < HU_DIFF_MIN)
        b2 = (cond & ~b1
              & (((h_upd[:, j] > h_upd_jb) & (tot_frz <= 0.0))
                 | (mflx_up[:, j] < MU_MIN)))
        jt0 = jnp.where(b1_low | b2, j + 1, jnp.where(b1, j, jt0))
        doit = doit & ~(b1 | b2)

    # zero the plume above the top; detrain what is left at the top
    for j in range(pver - 1, msg - 1, -1):
        m1 = (j >= lel) & (j <= jt0) & lm_pos
        mflx_up = mflx_up.at[:, j].set(
            jnp.where(m1, 0.0, mflx_up[:, j]))
        entr_up = entr_up.at[:, j].set(
            jnp.where(m1, 0.0, entr_up[:, j]))
        detr_up = detr_up.at[:, j].set(
            jnp.where(m1, 0.0, detr_up[:, j]))
        h_upd = h_upd.at[:, j].set(jnp.where(m1, h_env[:, j],
                                             h_upd[:, j]))
        m2 = (j == jt0) & lm_pos
        jp1 = min(j + 1, pver - 1)  # m2 empty at j=pver-1 (jt<jb<=pver-1)
        detr_up = detr_up.at[:, j].set(
            jnp.where(m2, mflx_up[:, jp1] / dz[:, j], detr_up[:, j]))
        entr_up = entr_up.at[:, j].set(
            jnp.where(m2, 0.0, entr_up[:, j]))
        mflx_up = mflx_up.at[:, j].set(
            jnp.where(m2, 0.0, mflx_up[:, j]))

    # LCL search with s/q recursion (bottom -> top, sequential)
    done = jnp.zeros(n, dtype=bool)
    for j in range(pver - 1, msg, -1):
        mjb = (j == jb) & lm_pos
        q_upd = q_upd.at[:, j].set(
            jnp.where(mjb, q_mid[:, j], q_upd[:, j]))
        s_upd = s_upd.at[:, j].set(
            jnp.where(mjb, (h_upd[:, j]
                            - zc["latvap"] * q_upd[:, j])
                      / zc["cpair"], s_upd[:, j]))
        m = (~done) & (j > jt0) & (j < jb) & lm_pos
        mf_s = jnp.where(mflx_up[:, j] != 0.0, mflx_up[:, j], 1.0)
        su_new = (mflx_up[:, j + 1] / mf_s * s_upd[:, j + 1]
                  + dz[:, j] / mf_s
                  * (entr_up[:, j] - detr_up[:, j]) * s_mid[:, j])
        qu_new = (mflx_up[:, j + 1] / mf_s * q_upd[:, j + 1]
                  + dz[:, j] / mf_s
                  * (entr_up[:, j] * q_mid[:, j]
                     - detr_up[:, j] * qst[:, j]))
        s_upd = s_upd.at[:, j].set(jnp.where(m, su_new, s_upd[:, j]))
        q_upd = q_upd.at[:, j].set(jnp.where(m, qu_new, q_upd[:, j]))
        tu = s_upd[:, j] - zc["grav"] / zc["cpair"] * z_int[:, j]
        _estu, qstu = qsat_hpa(
            tu, (p_mid[:, j] + p_mid[:, j - 1]) / 2.0)
        lclm = m & (q_upd[:, j] >= qstu)
        jlcl = jnp.where(lclm, j, jlcl)
        done = done | lclm

    # wet-adiabatic s/q between the top and the LCL
    m = ((ki > jt0[:, None]) & (ki <= jlcl[:, None])
         & lm_pos[:, None])
    s_upd = jnp.where(m, s_int + (h_upd - hsthat)
                      / (zc["cpair"] * (1.0 + gamhat)), s_upd)
    q_upd = jnp.where(m, qsthat + gamhat * (h_upd - hsthat)
                      / (zc["latvap"] * (1.0 + gamhat)), q_upd)

    # condensation rate in the updraft
    for j in range(pver - 1, msg, -1):
        m = lm_pos & (j >= jt0) & (j < jb)
        cu_new = (((mflx_up[:, j] * s_upd[:, j]
                    - mflx_up[:, j + 1] * s_upd[:, j + 1]) / dz[:, j]
                   - (entr_up[:, j] - detr_up[:, j]) * s_mid[:, j])
                  / (zc["latvap"] / zc["cpair"]))
        cu_new = jnp.where(j == jt0, 0.0, cu_new)
        cu = cu.at[:, j].set(
            jnp.where(m, jnp.maximum(0.0, cu_new), cu[:, j]))

    # in-plume liquid, rain production, total precip (no-microp path)
    rprd = jnp.zeros((n, pver))
    for j in range(pver - 1, msg, -1):
        m = ((j >= jt0) & (j < jb) & lm_pos
             & (mflx_up[:, j] >= 0.0))
        pos = m & (mflx_up[:, j] > 0.0)
        mf_s = jnp.where(mflx_up[:, j] != 0.0, mflx_up[:, j], 1.0)
        ql1 = (1.0 / mf_s
               * (mflx_up[:, j + 1] * ql[:, j + 1]
                  - dz[:, j] * detr_up[:, j] * ql[:, j + 1]
                  + dz[:, j] * cu[:, j]))
        ql = ql.at[:, j].set(
            jnp.where(pos, ql1 / (1.0 + dz[:, j] * c0mask),
                      jnp.where(m, 0.0, ql[:, j])))
        totpcp = totpcp + jnp.where(
            m, dz[:, j] * (cu[:, j] - detr_up[:, j] * ql[:, j + 1]),
            0.0)
        rprd = rprd.at[:, j].set(
            jnp.where(m, c0mask * mflx_up[:, j] * ql[:, j],
                      rprd[:, j]))

    (jt0, jd, mflx_dn, entr_dn, s_dnd, q_dnd, h_dnd, q_dnd_sat, evp,
     totevp) = zm_downdraft_properties(
        jb, jt0, j0, z_int, dz, s_mid, q_mid, h_env, lambda_max,
        qsthat, hsthat, gamhat, rprd, mflx_up, msg, zc, zp,
        s_dnd, q_dnd, q_dnd_sat, h_dnd)

    totpcp = jnp.maximum(totpcp, 0.0)
    totevp = jnp.maximum(totevp, 0.0)

    # downdraft strength consistent with precip availability
    posm = (totevp > 0.0) & (totpcp > 0.0)
    fac = jnp.minimum(
        1.0, totpcp / jnp.where(posm, totevp + totpcp, 1.0))
    m = (ki >= msg + 1)
    mflx_dn = jnp.where(m, jnp.where(posm[:, None],
                                     mflx_dn * fac[:, None], 0.0),
                        mflx_dn)
    entr_dn = jnp.where(m, jnp.where(posm[:, None],
                                     entr_dn * fac[:, None], 0.0),
                        entr_dn)
    evp = jnp.where(m, jnp.where(posm[:, None],
                                 evp * fac[:, None], 0.0), evp)
    rprd = rprd.at[:, msg + 1:].set(rprd[:, msg + 1:]
                                    - evp[:, msg + 1:])

    # net precipitation flux across interfaces
    pflx = jnp.zeros((n, pver + 1))
    for j in range(1, pver + 1):
        pflx = pflx.at[:, j].set(pflx[:, j - 1]
                                 + rprd[:, j - 1] * dz[:, j - 1])

    mflx_net = mflx_up + mflx_dn

    return dict(jt=jt0, jlcl=jlcl, j0=j0, jd=jd, mflx_up=mflx_up,
                entr_up=entr_up, detr_up=detr_up, mflx_dn=mflx_dn,
                entr_dn=entr_dn, mflx_net=mflx_net, s_upd=s_upd,
                q_upd=q_upd, ql=ql, s_dnd=s_dnd, q_dnd=q_dnd,
                qst=qst, cu=cu, evp=evp, pflx=pflx, rprd=rprd, dz=dz)


def zm_closure(lcl, lel, jt0, mx, dsubcld, z_int, p_mid, p_del, t_mid,
               s_mid, q_mid, qs, ql, s_int, q_int, t_pcl_lcl, t_pcl,
               q_pcl_sat, s_upd, q_upd, mflx_net, detr_up, mflx_up,
               mflx_dn, q_dnd, s_dnd, cape, cape_threshold, msg, zc,
               zp):
    """Z02 quasi-equilibrium closure: CAPE consumption rate per unit
    cloud-base mass flux -> cloud-base mass flux."""
    n, pver = t_mid.shape
    rows = jnp.arange(n)
    ki = jnp.arange(pver)[None, :]
    eps = zc["epsilo"]
    lat_cp = zc["latvap"] / zc["cpair"]
    cp_lat = zc["cpair"] / zc["latvap"]

    kmin0 = int(jnp.min(lel))
    kmax0 = int(jnp.max(mx)) - 1

    # sub-cloud tendencies
    p_mx = p_mid[rows, mx]
    q_mx = q_mid[rows, mx]
    t_mx = t_mid[rows, mx]
    eb = p_mx * q_mx / (eps + q_mx)
    dtbdt = (1.0 / dsubcld) * (
        mflx_up[rows, mx] * (s_int[rows, mx] - s_upd[rows, mx])
        + mflx_dn[rows, mx] * (s_int[rows, mx] - s_dnd[rows, mx]))
    dqbdt = (1.0 / dsubcld) * (
        mflx_up[rows, mx] * (q_int[rows, mx] - q_upd[rows, mx])
        + mflx_dn[rows, mx] * (q_int[rows, mx] - q_dnd[rows, mx]))
    debdt = eps * p_mx / (eps + q_mx)**2 * dqbdt
    dtldt = (-2840.0 * (3.5 / t_mx * dtbdt - debdt / eb)
             / (3.5 * jnp.log(t_mx) - jnp.log(eb) - 4.805)**2)

    # free-troposphere tendencies (vectorized over k)
    dtmdt = jnp.zeros((n, pver))
    dqmdt = jnp.zeros((n, pver))
    kv = slice(msg, pver - 1)
    kiv = ki[:, kv]
    mtop = kiv == jt0[:, None]
    dtmdt_top = (1.0 / p_del[:, kv]) * (
        mflx_up[:, msg + 1:] * (s_upd[:, msg + 1:]
                                - s_int[:, msg + 1:]
                                - lat_cp * ql[:, msg + 1:])
        + mflx_dn[:, msg + 1:] * (s_dnd[:, msg + 1:]
                                  - s_int[:, msg + 1:]))
    dqmdt_top = (1.0 / p_del[:, kv]) * (
        mflx_up[:, msg + 1:] * (q_upd[:, msg + 1:]
                                - q_int[:, msg + 1:]
                                + ql[:, msg + 1:])
        + mflx_dn[:, msg + 1:] * (q_dnd[:, msg + 1:]
                                  - q_int[:, msg + 1:]))
    mbelow = (kiv > jt0[:, None]) & (kiv < mx[:, None])
    qlb = BETA * ql[:, kv] + (1.0 - BETA) * ql[:, msg + 1:]
    dtmdt_bel = ((mflx_net[:, kv] * (s_int[:, kv] - s_mid[:, kv])
                  + mflx_net[:, msg + 1:] * (s_mid[:, kv]
                                             - s_int[:, msg + 1:]))
                 / p_del[:, kv]
                 - lat_cp * detr_up[:, kv] * qlb)

    def _qflux(mf, qq, ss, k_this):
        return mf * (qq - q_int[:, k_this]
                     + cp_lat * (ss - s_mid[:, kv]))

    dqmdt_bel = ((_qflux(mflx_up[:, msg + 1:], q_upd[:, msg + 1:],
                         s_upd[:, msg + 1:], slice(msg + 1, pver))
                  - _qflux(mflx_up[:, kv], q_upd[:, kv],
                           s_upd[:, kv], kv)
                  + _qflux(mflx_dn[:, msg + 1:], q_dnd[:, msg + 1:],
                           s_dnd[:, msg + 1:], slice(msg + 1, pver))
                  - _qflux(mflx_dn[:, kv], q_dnd[:, kv],
                           s_dnd[:, kv], kv)) / p_del[:, kv]
                 + detr_up[:, kv] * qlb)
    dtmdt = dtmdt.at[:, kv].set(
        jnp.where(mtop, dtmdt_top,
                  jnp.where(mbelow, dtmdt_bel, 0.0)))
    dqmdt = dqmdt.at[:, kv].set(
        jnp.where(mtop, dqmdt_top,
                  jnp.where(mbelow, dqmdt_bel, 0.0)))

    # dboydt (integrand of the CAPE change) - vectorized over k
    kappa = zc["rdair"] / zc["cpair"]
    mdry = (ki > lcl[:, None]) & (ki < mx[:, None]) & (ki >= msg)
    mwet = (ki >= lel[:, None]) & (ki <= lcl[:, None]) & (ki >= msg)
    pfac = (1000.0 / p_mid)**kappa
    thetavm = t_mid * pfac * (1.0 + 0.608 * q_mid)
    q_mx_c = q_mx[:, None]
    thetavp_dry = t_pcl * pfac * (1.0 + 0.608 * q_mx_c)
    denvm = (dtmdt / t_mid
             + 0.608 / (1.0 + 0.608 * q_mid) * dqmdt)
    dbdt_b = (dtbdt / t_mx)[:, None] \
        + (0.608 / (1.0 + 0.608 * q_mx) * dqbdt)[:, None]
    dboydt_dry = ((dbdt_b - denvm) * zc["grav"]
                  * thetavp_dry / thetavm)

    t_pcl_s = jnp.where(t_pcl != 0.0, t_pcl, 1.0)  # masked-off lanes
    thetavp_wet = (t_pcl * pfac
                   * (1.0 + 1.608 * q_pcl_sat - q_mx_c))
    dqsdtp = (q_pcl_sat * (1.0 + q_pcl_sat / eps) * eps
              * zc["latvap"] / (zc["rdair"] * t_pcl_s**2))
    dtpdt = (t_pcl / (1.0 + lat_cp
                      * (dqsdtp - q_pcl_sat / t_pcl_s))
             * ((dtbdt / t_mx)[:, None]
                + lat_cp * ((dqbdt / t_pcl_lcl)[:, None]
                            - (q_mx / t_pcl_lcl**2
                               * dtldt)[:, None])))
    dboydt_wet = ((dtpdt / t_pcl_s
                   + 1.0 / (1.0 + 1.608 * q_pcl_sat - q_mx_c)
                   * (1.608 * dqsdtp * dtpdt - dqbdt[:, None]))
                  - denvm) * zc["grav"] * thetavp_wet / thetavm
    dboydt = jnp.where(mdry, dboydt_dry,
                       jnp.where(mwet, dboydt_wet, 0.0))

    # vertical integral of the buoyancy change (Fortran loop order)
    dadt = jnp.zeros(n)
    for j in range(kmin0, kmax0 + 1):
        m = (j >= lel) & (j <= mx - 1)
        dadt = dadt + jnp.where(
            m, dboydt[:, j] * (z_int[:, j] - z_int[:, j + 1]), 0.0)

    dltaa = -1.0 * (cape - cape_threshold)
    cbmf = jnp.where(dadt != 0.0,
                     jnp.maximum(dltaa / zp["tau"]
                                 / jnp.where(dadt != 0.0, dadt, 1.0),
                                 0.0), 0.0)
    return cbmf


def zm_calc_output_tend(jt0, mx, dsubcld, p_del, s_int, q_int, s_upd,
                        q_upd, mflx_up, detr_up, mflx_dn, s_dnd,
                        q_dnd, ql, evp, cu, msg, zc):
    """Final dsdt/dqdt/dl tendencies (zm_microp=False branch)."""
    n, pver = p_del.shape
    ki = jnp.arange(pver)[None, :]
    lat_cp = zc["latvap"] / zc["cpair"]
    dsdt = jnp.zeros((n, pver))
    dqdt = jnp.zeros((n, pver))
    dl = jnp.zeros((n, pver))

    ktm0 = int(jnp.min(jt0))
    kbm0 = int(jnp.min(mx))

    # free troposphere (vectorized; unconditional across gathered
    # columns exactly as the Fortran k=ktm..pver-1 / i loops)
    kv = slice(ktm0, pver - 1)
    kp = slice(ktm0 + 1, pver)
    emc = -cu[:, kv] + evp[:, kv]
    dsdt = dsdt.at[:, kv].set(
        -lat_cp * emc
        + (mflx_up[:, kp] * (s_upd[:, kp] - s_int[:, kp])
           - mflx_up[:, kv] * (s_upd[:, kv] - s_int[:, kv])
           + mflx_dn[:, kp] * (s_dnd[:, kp] - s_int[:, kp])
           - mflx_dn[:, kv] * (s_dnd[:, kv] - s_int[:, kv]))
        / p_del[:, kv])
    dqdt = dqdt.at[:, kv].set(
        emc
        + (mflx_up[:, kp] * (q_upd[:, kp] - q_int[:, kp])
           - mflx_up[:, kv] * (q_upd[:, kv] - q_int[:, kv])
           + mflx_dn[:, kp] * (q_dnd[:, kp] - q_int[:, kp])
           - mflx_dn[:, kv] * (q_dnd[:, kv] - q_int[:, kv]))
        / p_del[:, kv])
    dl = dl.at[:, kv].set(detr_up[:, kv] * ql[:, kp])

    # at and below cloud base (sequential: k > mx copies k-1)
    for j in range(kbm0, pver):
        mmx = j == mx
        dsdt = dsdt.at[:, j].set(jnp.where(
            mmx, (1.0 / dsubcld)
            * (-mflx_up[:, j] * (s_upd[:, j] - s_int[:, j])
               - mflx_dn[:, j] * (s_dnd[:, j] - s_int[:, j])),
            dsdt[:, j]))
        dqdt = dqdt.at[:, j].set(jnp.where(
            mmx, (1.0 / dsubcld)
            * (-mflx_up[:, j] * (q_upd[:, j] - q_int[:, j])
               - mflx_dn[:, j] * (q_dnd[:, j] - q_int[:, j])),
            dqdt[:, j]))
        mgt = j > mx
        dsdt = dsdt.at[:, j].set(jnp.where(mgt, dsdt[:, j - 1],
                                           dsdt[:, j]))
        dqdt = dqdt.at[:, j].set(jnp.where(mgt, dqdt[:, j - 1],
                                           dqdt[:, j]))
    return dsdt, dqdt, dl


def zm_conv_main(t_mid, q_mid_in, omega, p_mid_in, p_int_in, p_del_in,
                 geos, z_mid_in, z_int_in, pbl_hgt, tpert, landfrac,
                 t_star, q_star, time_step, is_first_step, zc, zp):
    """ZM deep convection main routine (zm_microp=False). Inputs are
    (ncol, pver)/(ncol, pver+1) C-ordered arrays, level 0 = top;
    pressures in Pa, z relative to the surface [m], geos [m2/s2].
    Returns a dict of all zm_conv_main outputs (0-based indices;
    gathered/scattered layout as in PORT_NOTES)."""
    t_mid = jnp.asarray(t_mid, dtype=jnp.float64)
    q_mid_in = jnp.asarray(q_mid_in, dtype=jnp.float64)
    omega = jnp.asarray(omega, dtype=jnp.float64)
    p_del_in = jnp.asarray(p_del_in, dtype=jnp.float64)
    ncol, pver = t_mid.shape
    msg = zp["limcnv"] - 1

    # local pressure [mb] and height [m] incl. surface elevation
    z_srf = jnp.asarray(geos) * (1.0 / zc["grav"])
    p_mid = jnp.asarray(p_mid_in) * 0.01
    p_int = jnp.asarray(p_int_in) * 0.01
    z_mid = jnp.asarray(z_mid_in) + z_srf[:, None]
    z_int = jnp.asarray(z_int_in) + z_srf[:, None]

    # PBL top index (descending k; the highest matching level wins)
    pbl_top = jnp.full(ncol, pver - 1, dtype=jnp.int64)
    for j in range(pver - 2, msg - 1, -1):
        m = (jnp.abs(z_mid[:, j] - z_srf - jnp.asarray(pbl_hgt))
             < (z_int[:, j] - z_int[:, j + 1]) * 0.5)
        pbl_top = jnp.where(m, j, pbl_top)

    s_mid = t_mid + (zc["grav"] / zc["cpair"]) * z_mid

    # dilute CAPE, call 1 (current state)
    r1 = compute_dilute_cape(q_mid_in, t_mid, z_mid, p_mid, p_int,
                             pbl_top, tpert, msg, zc, zp,
                             calc_msemax_klev=True)
    cape = r1["cape"]
    msemax_klev = r1["msemax_klev"]

    # DCAPE trigger: call 2 on the previous state, frozen launch level
    dcape = jnp.zeros(ncol)
    if (not is_first_step) and zp["trig_dcape"]:
        r2 = compute_dilute_cape(q_star, t_star, z_mid, p_mid, p_int,
                                 pbl_top, tpert, msg, zc, zp,
                                 calc_msemax_klev=False,
                                 prev_msemax_klev=msemax_klev)
        dcape = (cape - r2["cape"]) / time_step

    # gather trigger (zm_get_gather_index)
    if zp["trig_dcape"] and not is_first_step:
        cape_threshold_loc = CAPE_THRESHOLD_NEW
        trig = (np.asarray(cape) > cape_threshold_loc) \
            & (np.asarray(dcape) > DCAPE_THRESHOLD)
    else:
        cape_threshold_loc = CAPE_THRESHOLD_OLD
        trig = np.asarray(cape) > cape_threshold_loc
    gidx = np.nonzero(trig)[0]
    lengath = int(gidx.size)

    # outputs with the exact Fortran no-op fill values
    out = dict(
        lengath=lengath,
        gather_index=np.full(ncol, -1, dtype=np.int64),
        msemax_klev_g=jnp.zeros(ncol, dtype=jnp.int64),
        jctop=jnp.full(ncol, pver - 1, dtype=jnp.int64),
        jcbot=jnp.zeros(ncol, dtype=jnp.int64),
        jt=jnp.full(ncol, pver - 1, dtype=jnp.int64),
        prec=jnp.zeros(ncol), heat=jnp.zeros((ncol, pver)),
        qtnd=jnp.zeros((ncol, pver)), cape=cape, dcape=dcape,
        mcon=jnp.zeros((ncol, pver + 1)),
        pflx=jnp.zeros((ncol, pver + 1)),
        zdu=jnp.zeros((ncol, pver)),
        mflx_up=jnp.zeros((ncol, pver)),
        entr_up=jnp.zeros((ncol, pver)),
        detr_up=jnp.zeros((ncol, pver)),
        mflx_dn=jnp.zeros((ncol, pver)),
        entr_dn=jnp.zeros((ncol, pver)),
        p_del=jnp.zeros((ncol, pver)), dsubcld=jnp.zeros(ncol),
        ql=jnp.zeros((ncol, pver)), rliq=jnp.zeros(ncol),
        rprd=jnp.zeros((ncol, pver)), dlf=jnp.zeros((ncol, pver)),
        cld_base_mass_flux=jnp.zeros(ncol))
    out["gather_index"][:lengath] = gidx
    if lengath == 0:
        return out

    # ---- gathered arrays ----
    g = gidx
    p_del_g = 0.01 * p_del_in[g]
    q_mid_g = q_mid_in[g]
    t_mid_g = t_mid[g]
    p_mid_g = p_mid[g]
    z_mid_g = z_mid[g]
    s_mid_g = s_mid[g]
    t_pcl_g = r1["parcel_temp"][g]
    z_int_g = z_int[g]
    q_pcl_sat_g = r1["parcel_qsat"][g]
    omega_g = omega[g]
    cape_g = cape[g]
    lcl_g = r1["lcl_klev"][g]
    lel_g = r1["eql_klev"][g]
    mx_g = msemax_klev[g]
    t_pcl_lcl_g = r1["lcl_temperature"][g]
    landfrac_g = jnp.asarray(landfrac)[g]
    pbl_top_g = pbl_top[g]
    tpert_g = jnp.asarray(tpert)[g]
    rows = jnp.arange(lengath)
    ki = jnp.arange(pver)[None, :]

    # sub-cloud layer pressure thickness
    dsubcld = jnp.sum(jnp.where((ki >= mx_g[:, None]) & (ki >= msg),
                                p_del_g, 0.0), axis=1)

    # interfacial (q,s) values
    s_int_g = _interface_interp(s_mid_g, msg)
    q_int_g = _interface_interp(q_mid_g, msg)

    cp = zm_cloud_properties(p_mid_g, z_mid_g, z_int_g, t_mid_g,
                             s_mid_g, s_int_g, q_mid_g, landfrac_g,
                             tpert_g, mx_g, lel_g, msg, zc, zp)
    jt_g = cp["jt"]
    dz = cp["dz"]
    mflx_up = cp["mflx_up"]
    entr_up = cp["entr_up"]
    detr_up = cp["detr_up"]
    mflx_dn = cp["mflx_dn"]
    entr_dn = cp["entr_dn"]
    mflx_net = cp["mflx_net"]
    ql_g = cp["ql"]
    cu_g = cp["cu"]
    evp_g = cp["evp"]
    pflx_g = cp["pflx"]
    rprd_g = cp["rprd"]

    # convert "per length" [1/m] -> "per pressure" [1/mb]
    conv = jnp.where(ki >= msg, dz / p_del_g, 0.0)
    keep = ki < msg
    detr_up = jnp.where(keep, detr_up, detr_up * conv)
    entr_up = jnp.where(keep, entr_up, entr_up * conv)
    entr_dn = jnp.where(keep, entr_dn, entr_dn * conv)
    cu_g = jnp.where(keep, cu_g, cu_g * conv)
    rprd_g = jnp.where(keep, rprd_g, rprd_g * conv)
    evp_g = jnp.where(keep, evp_g, evp_g * conv)

    # CAPE closure -> cloud-base mass flux
    cbmf = zm_closure(lcl_g, lel_g, jt_g, mx_g, dsubcld, z_int_g,
                      p_mid_g, p_del_g, t_mid_g, s_mid_g, q_mid_g,
                      cp["qst"], ql_g, s_int_g, q_int_g, t_pcl_lcl_g,
                      t_pcl_g, q_pcl_sat_g, cp["s_upd"], cp["q_upd"],
                      mflx_net, detr_up, mflx_up, mflx_dn,
                      cp["q_dnd"], cp["s_dnd"], cape_g,
                      cape_threshold_loc, msg, zc, zp)

    # limit to the theoretical upper bound
    mflx_up_max = jnp.max(
        jnp.where(ki >= msg + 1, mflx_up / p_del_g, 0.0), axis=1)
    mflx_up_max = jnp.maximum(mflx_up_max, 0.0)
    cbmf = jnp.where(
        mflx_up_max > 0.0,
        jnp.minimum(cbmf, 1.0 / (time_step
                                 * jnp.where(mflx_up_max > 0.0,
                                             mflx_up_max, 1.0))),
        0.0)
    if zp["clos_dyn_adj"]:
        cbmf = jnp.maximum(
            cbmf - omega_g[rows, pbl_top_g] * 0.01, 0.0)

    # no deep convection within the PBL
    if zp["no_deep_pbl"]:
        zjt = jnp.asarray(z_mid_in)[g][rows, jt_g]
        cbmf = jnp.where(zjt < jnp.asarray(pbl_hgt)[g], 0.0, cbmf)

    # scale by the cloud-base mass flux
    sc = cbmf[:, None]
    m = ki >= msg
    mflx_up = jnp.where(m, mflx_up * sc, mflx_up)
    mflx_dn = jnp.where(m, mflx_dn * sc, mflx_dn)
    mflx_net = jnp.where(m, mflx_net * sc, mflx_net)
    detr_up = jnp.where(m, detr_up * sc, detr_up)
    entr_up = jnp.where(m, entr_up * sc, entr_up)
    entr_dn = jnp.where(m, entr_dn * sc, entr_dn)
    rprd_g = jnp.where(m, rprd_g * sc, rprd_g)
    cu_g = jnp.where(m, cu_g * sc, cu_g)
    evp_g = jnp.where(m, evp_g * sc, evp_g)
    kint = jnp.arange(pver + 1)[None, :]
    pflx_g = jnp.where((kint >= msg + 1) & (kint <= pver),
                       pflx_g * sc * 100.0 / zc["grav"], pflx_g)

    # output tendencies
    dsdt, dqdt, dl_g = zm_calc_output_tend(
        jt_g, mx_g, dsubcld, p_del_g, s_int_g, q_int_g, cp["s_upd"],
        cp["q_upd"], mflx_up, detr_up, mflx_dn, cp["s_dnd"],
        cp["q_dnd"], ql_g, evp_g, cu_g, msg, zc)

    # ---- scatter ----
    gr = jnp.asarray(g)
    mk = jnp.arange(pver)[None, :] >= msg
    q_mid_new = q_mid_in[g] + time_step * dqdt
    qtnd = out["qtnd"].at[gr].set(jnp.where(mk, dqdt, 0.0))
    heat = out["heat"].at[gr].set(jnp.where(mk, dsdt * zc["cpair"],
                                            0.0))
    rprd = out["rprd"].at[gr].set(jnp.where(mk, rprd_g, 0.0))
    zdu = out["zdu"].at[gr].set(jnp.where(mk, detr_up, 0.0))
    mcon = out["mcon"].at[gr, :pver].set(jnp.where(mk, mflx_net, 0.0))
    dlf = out["dlf"].at[gr].set(jnp.where(mk, dl_g, 0.0))
    pflx = out["pflx"].at[gr, :pver].set(
        jnp.where(mk, pflx_g[:, :pver], 0.0))
    pflx = pflx.at[gr, pver].set(pflx_g[:, pver])
    ql_out = out["ql"].at[gr].set(jnp.where(mk, ql_g, 0.0))
    jctop = out["jctop"].at[gr].set(jt_g)
    jcbot = out["jcbot"].at[gr].set(mx_g)

    # precip from the change in water vapor minus detrained liquid
    # (all columns; descending-k accumulation order preserved)
    dq_full = jnp.zeros((ncol, pver)).at[gr].set(
        jnp.where(mk, q_mid_new - q_mid_in[g], 0.0))
    prec = jnp.zeros(ncol)
    for j in range(pver - 1, msg - 1, -1):
        prec = prec - p_del_in[:, j] * dq_full[:, j] \
            - p_del_in[:, j] * dlf[:, j] * time_step
    prec = (1.0 / zc["grav"]) * jnp.maximum(prec, 0.0) \
        / time_step / 1000.0

    # reserved liquid (ascending-k accumulation order preserved)
    rliq = jnp.zeros(ncol)
    for j in range(pver):
        rliq = rliq + dlf[:, j] * p_del_in[:, j] / zc["grav"]
    rliq = rliq / 1000.0

    # gathered-layout outputs (rows >= lengath keep the fill values)
    out["msemax_klev_g"] = out["msemax_klev_g"].at[rows].set(mx_g)
    out["jt"] = out["jt"].at[rows].set(jt_g)
    out["mflx_up"] = out["mflx_up"].at[rows].set(mflx_up)
    out["entr_up"] = out["entr_up"].at[rows].set(entr_up)
    out["detr_up"] = out["detr_up"].at[rows].set(detr_up)
    out["mflx_dn"] = out["mflx_dn"].at[rows].set(mflx_dn)
    out["entr_dn"] = out["entr_dn"].at[rows].set(entr_dn)
    out["p_del"] = out["p_del"].at[rows].set(p_del_g)
    out["dsubcld"] = out["dsubcld"].at[rows].set(dsubcld)
    out["ql"] = ql_out
    out["cld_base_mass_flux"] = \
        out["cld_base_mass_flux"].at[rows].set(cbmf)

    out.update(jctop=jctop, jcbot=jcbot, prec=prec, heat=heat,
               qtnd=qtnd, mcon=mcon, pflx=pflx, zdu=zdu, rliq=rliq,
               rprd=rprd, dlf=dlf)
    return out


# ---------------------------------------------------------------------------
# Precipitation evaporation / snow production (zm_conv_evap)
# ---------------------------------------------------------------------------
def cldfrc_fice(t_mid, zc):
    """cloud_fraction.F90 cldfrc_fice: fraction of condensate in the
    ice phase and the snow fraction for convection, both piecewise
    linear in temperature (top_lev = 1: all levels computed; see
    PORT_NOTES). Returns (fice, fsnow)."""
    t_mid = jnp.asarray(t_mid, dtype=jnp.float64)
    tmax_fice = zc["tfreez"] - 10.0
    tmin_fice = tmax_fice - 30.0
    tmax_fsnow = zc["tfreez"]
    tmin_fsnow = zc["tfreez"] - 5.0
    fice = jnp.where(
        t_mid > tmax_fice, 0.0,
        jnp.where(t_mid < tmin_fice, 1.0,
                  (tmax_fice - t_mid) / (tmax_fice - tmin_fice)))
    fsnow = jnp.where(
        t_mid > tmax_fsnow, 0.0,
        jnp.where(t_mid < tmin_fsnow, 1.0,
                  (tmax_fsnow - t_mid) / (tmax_fsnow - tmin_fsnow)))
    return fice, fsnow


def zm_conv_evap(p_mid, p_del, t_mid, q_mid, prdprec, cldfrc, prec_in,
                 time_step, zc, zp, prdsnow=None):
    """zm_conv.F90 zm_conv_evap: below-cloud Sundqvist-type
    evaporation of convective precip, snow melt/production with the
    latent heat of fusion, and the surface precip/snow rates.

    p_mid/p_del [Pa], t_mid [K], q_mid [kg/kg], prdprec [kg/kg/s]
    (rain production, rprd), cldfrc (cloud fraction), prec_in [m/s]
    (the intent(inout) prec from zm_conv_main), time_step [s].
    zp needs ke, old_snow (zm_microp is False in this scope; prdsnow
    defaults to zero, exactly the .not. zm_microp branch).
    tend_s/tend_q are pure outputs (see PORT_NOTES). Returns a dict:
    tend_s, tend_q, tend_s_snwprd, tend_s_snwevmlt [J/kg/s | kg/kg/s],
    prec, snow [m/s], ntprprd, ntsnprd [kg/kg/s], flxprec, flxsnow
    [kg/m2/s] (ncol, pver+1)."""
    p_mid = jnp.asarray(p_mid, dtype=jnp.float64)
    p_del = jnp.asarray(p_del, dtype=jnp.float64)
    t_mid = jnp.asarray(t_mid, dtype=jnp.float64)
    q_mid = jnp.asarray(q_mid, dtype=jnp.float64)
    prdprec = jnp.asarray(prdprec, dtype=jnp.float64)
    cldfrc = jnp.asarray(cldfrc, dtype=jnp.float64)
    ncol, pver = t_mid.shape
    grav = zc["grav"]
    latice = zc["latice"]
    latvap = zc["latvap"]
    cpair = zc["cpair"]
    tfreez = zc["tfreez"]
    old_snow = bool(zp["old_snow"])
    ke = zp["ke"]
    pergro_active = False

    if prdsnow is None:
        prdsnow = jnp.zeros((ncol, pver))   # .not. zm_microp
    else:
        prdsnow = jnp.asarray(prdsnow, dtype=jnp.float64)

    # convert input precip to kg/m2/s
    prec = jnp.asarray(prec_in, dtype=jnp.float64) * 1000.0

    # saturation (mixed-phase table qsat, as wv_saturation's qsat)
    qs = _qsat_pa(t_mid, p_mid)["qs"]

    # ice/snow fraction in rain production
    _fice, fsnow_conv = cldfrc_fice(t_mid, zc)

    flxprec = jnp.zeros((ncol, pver + 1))
    flxsnow = jnp.zeros((ncol, pver + 1))
    evpvint = jnp.zeros(ncol)
    tend_s = jnp.zeros((ncol, pver))
    tend_q = jnp.zeros((ncol, pver))
    tend_s_snwprd = jnp.zeros((ncol, pver))
    tend_s_snwevmlt = jnp.zeros((ncol, pver))
    ntprprd = jnp.zeros((ncol, pver))
    ntsnprd = jnp.zeros((ncol, pver))

    for k in range(pver):
        tk = t_mid[:, k]
        dpk = p_del[:, k]
        fpk = flxprec[:, k]
        fsk = flxsnow[:, k]
        warm = tk > tfreez

        # melt snow falling into the layer
        if old_snow:
            flxsntm = jnp.where(warm, 0.0, fsk)
            snowmlt = jnp.where(warm, fsk * grav / dpk, 0.0)
        else:
            # make sure melting snow doesn't cool below the threshold
            dum0 = -latice / cpair * fsk * grav / dpk * time_step
            full_melt_too_cold = (tk + dum0) <= tfreez
            denom = jnp.where(full_melt_too_cold & warm,
                              fsk * grav / dpk, 1.0)
            dum_lim = (tk - tfreez) * cpair / latice / time_step / denom
            dum_lim = jnp.clip(dum_lim, 0.0, 1.0)
            dum = jnp.where(full_melt_too_cold, dum_lim, 1.0) * OMSM
            flxsntm = jnp.where(warm, fsk * (1.0 - dum), fsk)
            snowmlt = jnp.where(warm, dum * fsk * grav / dpk, 0.0)

        # relative humidity depression must be > 0 for evaporation
        evplimit = jnp.maximum(1.0 - q_mid[:, k] / qs[:, k], 0.0)
        evpprec = (ke * (1.0 - cldfrc[:, k]) * evplimit
                   * jnp.sqrt(fpk))

        # don't supersaturate; don't evaporate more than falls in;
        # don't exceed the remaining input precipitation
        evplimit = jnp.maximum(0.0, (qs[:, k] - q_mid[:, k]) / time_step)
        evplimit = jnp.minimum(evplimit, fpk * grav / dpk)
        evplimit = jnp.minimum(evplimit, (prec - evpvint) * grav / dpk)
        evpprec = jnp.minimum(evplimit, evpprec)
        if not old_snow:
            evpprec = jnp.maximum(0.0, evpprec) * OMSM

        # snow evaporation from the post-melt snow fraction
        fp_pos = fpk > 0.0
        work1 = jnp.clip(flxsntm / jnp.where(fp_pos, fpk, 1.0), 0.0, 1.0)
        if not old_snow:
            work1 = jnp.where(prdsnow[:, k] > prdprec[:, k], 1.0, work1)
        evpsnow = jnp.where(fp_pos, evpprec * work1, 0.0)

        # vertically integrated evaporation
        evpvint = evpvint + evpprec * dpk / grav

        # net precip production
        ntp = prdprec[:, k] - evpprec
        ntprprd = ntprprd.at[:, k].set(ntp)

        # net snow production
        if old_snow:
            if pergro_active:
                work1b = jnp.clip(fsk / (fpk + 8.64e-11), 0.0, 1.0)
            else:
                work1b = jnp.where(fp_pos,
                                   jnp.clip(fsk / jnp.where(fp_pos, fpk,
                                                            1.0),
                                            0.0, 1.0), 0.0)
            work2 = jnp.maximum(fsnow_conv[:, k], work1b)
            work2 = jnp.where(snowmlt > 0.0, 0.0, work2)
            nts = prdprec[:, k] * work2 - evpsnow - snowmlt
            tssp = prdprec[:, k] * work2 * latice
            tsse = -(evpsnow + snowmlt) * latice
        else:
            sink = jnp.minimum(fsk * grav / dpk, evpsnow + snowmlt)
            nts = prdsnow[:, k] - sink
            tssp = prdsnow[:, k] * latice
            tsse = -sink * latice
        ntsnprd = ntsnprd.at[:, k].set(nts)
        tend_s_snwprd = tend_s_snwprd.at[:, k].set(tssp)
        tend_s_snwevmlt = tend_s_snwevmlt.at[:, k].set(tsse)

        # precipitation fluxes, protected against rounding error
        flxprec = flxprec.at[:, k + 1].set(
            jnp.maximum(fpk + ntp * dpk / grav, 0.0))
        flxsnow = flxsnow.at[:, k + 1].set(
            jnp.maximum(fsk + nts * dpk / grav, 0.0))

        # heating/cooling and moistening due to evaporation
        if old_snow:
            tend_s = tend_s.at[:, k].set(-evpprec * latvap
                                         + nts * latice)
        else:
            tend_s = tend_s.at[:, k].set(-evpprec * latvap + tsse)
        tend_q = tend_q.at[:, k].set(evpprec)

    # protect against rounding error: enforce flxsnow <= flxprec at
    # the surface (only possible when prdsnow is active, i.e. with
    # zm_microphysics; unreachable in this zm_microp=False scope)
    if not old_snow:
        excess = flxsnow[:, pver] > flxprec[:, pver]
        dum = jnp.where(excess,
                        (flxsnow[:, pver] - flxprec[:, pver]) * grav,
                        0.0)
        for k in range(pver - 1, -1, -1):
            m = (ntsnprd[:, k] > ntprprd[:, k]) & (dum > 0.0)
            adj = dum / p_del[:, k]
            ntsnprd = ntsnprd.at[:, k].set(
                jnp.where(m, ntsnprd[:, k] - adj, ntsnprd[:, k]))
            tend_s_snwevmlt = tend_s_snwevmlt.at[:, k].set(
                jnp.where(m, tend_s_snwevmlt[:, k] - adj * latice,
                          tend_s_snwevmlt[:, k]))
            tend_s = tend_s.at[:, k].set(
                jnp.where(m, tend_s[:, k] - adj * latice, tend_s[:, k]))
            dum = jnp.where(m, 0.0, dum)
        flxsnow = flxsnow.at[:, pver].set(
            jnp.where(excess, flxprec[:, pver], flxsnow[:, pver]))

    # output precipitation rates [m/s]
    prec_out = flxprec[:, pver] / 1000.0
    snow_out = flxsnow[:, pver] / 1000.0

    return dict(tend_s=tend_s, tend_q=tend_q,
                tend_s_snwprd=tend_s_snwprd,
                tend_s_snwevmlt=tend_s_snwevmlt,
                prec=prec_out, snow=snow_out,
                ntprprd=ntprprd, ntsnprd=ntsnprd,
                flxprec=flxprec, flxsnow=flxsnow)
