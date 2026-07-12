"""Port of the ZM deep-convection dilute CAPE core,
eam/src/physics/cam/zm/zm_conv_cape.F90 (compute_dilute_cape =
find_mse_max + compute_dilute_parcel + compute_cape_from_parcel) with
its zm_conv_util.F90 helpers (entropy, ientropy, qsat_hPa).

PORT_NOTES
----------
- Index convention: everything here is 0-based, level 0 = model top,
  level pver-1 = surface. Fortran midpoint k maps to python k-1, so
  the integer level inputs/outputs (pblt, prev_msemax_klev,
  msemax_klev, lcl_klev, eql_klev) are Fortran values minus 1.
  num_msg is unchanged: it is the *count* of excluded top levels
  (Fortran loops k = pver..num_msg+1 == python j = pver-1..num_msg).
- Pressures are in hPa (mb), exactly as zm_conv.F90 passes them
  (p_mid_in * 0.01); zmid is height above sea level [m] (z + z_srf).
- Constants (make_zm_const) mirror zm_conv_types
  zm_const_set_to_global(): every value derives from
  share/util/shr_const_mod.F90 via the physconst derivations, EXCEPT
  zvir = 1.608 exactly, hardcoded there "to avoid non-BFB diffs"
  (NOT rh2o/rdair - 1).
- Parameters (make_zm_param) default to EAMv3 phys="default" values
  from bld/namelist_files/namelist_defaults_eam.xml:
  dmpdz = -0.7e-3 [1/m] (zmconv_dmpdz), tiedke_add = 0.8
  (zmconv_tiedke_add), tpert_fac = 2.0 (zmconv_tp_fac),
  mx_bot_lyr_adj = 1, num_cin = 1 (zmconv_cape_cin),
  tpert_fix / trig_ull / trig_dcape = True.
- ientropy is the Numerical-Recipes Brent inversion of entropy(T);
  ported with the *exact* Fortran update order (re-bracket, best-root
  swap, mid-body convergence exit, inverse-quadratic/secant/bisection
  choice, sign(tolerance, xm) minimum step) as a masked iteration:
  converged lanes freeze, so each lane follows the Fortran trajectory
  point by point. tol = 2*3e-8*|b| + 0.5*1e-3 (~5e-4 K). LOOPMAX =
  100; non-convergence raises (the Fortran endrun's). Masked-off lanes
  are fed self-consistent environment values so they converge quickly
  and never poison the run.
- Saturation comes from wv_saturation's qsat_water via the qsat_hPa
  hPa<->Pa wrapper, i.e. the already-ported eam_jax.wv_sat.qsat_water
  (GoffGratch default scheme, the wv_sat_init default).
- The entrainment loop, the LCL interpolation trigger (first level
  where qsmix <= qtmix on ascent), the two-pass (nit_lheat = 2)
  condensation/freezing adjustment with rainout above lwmax = 1e-3,
  the one-off/continual freezing branches at tfreez, and the
  num_cin-deep equilibrium-level bookkeeping all replicate the exact
  Fortran loop and branch order; sequential level sweeps are python
  loops with per-column masks, and the CAPE integral accumulates in
  the same ascending-k order (no jnp.sum, which reorders).
- q_mx/t_mx follow the harness driver's "always present" semantics:
  they are returned (saved at the launch level) on every call where
  use_input_tq_mx is False. This only affects those two outputs.
- lcl_pressure_threshold = ull_upper_launch_pressure = 600 hPa and
  pergro_rhd_threshold = -1e-4 are zm_conv_cape module parameters.
  pergro_active mirrors the PERGRO ifdef (goldens use False).
- Everything is float64, eager-mode (the Brent early-exit and level
  skipping use concrete booleans; not jit-compatible as written).
"""

import jax.numpy as jnp

from .constants import (CPAIR, CPWV, EPSILO, GRAVIT, LATICE, LATVAP,
                        RAIR, RH2O, SHR_CONST_CPFW, SHR_CONST_PI,
                        TMELT)
from .wv_sat import qsat_water

# zm_conv_cape module parameters
LCL_PRESSURE_THRESHOLD = 600.0     # [hPa] LCL above this -> no CAPE
ULL_UPPER_LAUNCH_PRESSURE = 600.0  # [hPa] ULL search upper limit
PERGRO_RHD_THRESHOLD = -1.0e-4
NIT_LHEAT = 2                      # condensation/freezing iterations
LWMAX = 1.0e-3                     # max retained condensate [kg/kg]

# zm_conv_util ientropy parameters
LOOPMAX = 100
TOL_COEFF = 0.001
TOL_EPS = 3.0e-8
PREF = 1000.0                      # [hPa] entropy reference pressure


def make_zm_const():
    """zm_const_t as filled by zm_const_set_to_global (physconst /
    shr_const_mod values; zvir hardcoded 1.608 there)."""
    return dict(pi=SHR_CONST_PI, grav=GRAVIT, rgrav=1.0 / GRAVIT,
                rdair=RAIR, rh2o=RH2O, zvir=1.608, cpair=CPAIR,
                cpwv=CPWV, cpliq=SHR_CONST_CPFW, tfreez=TMELT,
                latvap=LATVAP, latice=LATICE, epsilo=EPSILO)


def make_zm_param(dmpdz=-0.7e-3, tiedke_add=0.8, tpert_fac=2.0,
                  mx_bot_lyr_adj=1, num_cin=1, tpert_fix=True,
                  trig_ull=True, trig_dcape=True):
    """CAPE-relevant subset of zm_param_t; defaults are the EAMv3
    phys="default" namelist values (see PORT_NOTES)."""
    return dict(dmpdz=dmpdz, tiedke_add=tiedke_add,
                tpert_fac=tpert_fac, mx_bot_lyr_adj=int(mx_bot_lyr_adj),
                num_cin=int(num_cin), tpert_fix=bool(tpert_fix),
                trig_ull=bool(trig_ull), trig_dcape=bool(trig_dcape))


def _sign(a, b):
    """Fortran sign(a, b): |a| with the sign of b (+ for b == 0)."""
    return jnp.where(b >= 0.0, jnp.abs(a), -jnp.abs(a))


# ---------------------------------------------------------------------------
# zm_conv_util helpers
# ---------------------------------------------------------------------------
def qsat_hpa(t, p):
    """qsat_hPa: wv_saturation qsat_water with hPa<->Pa translation.
    Returns (es [hPa], qm [kg/kg])."""
    r = qsat_water(t, jnp.asarray(p, dtype=jnp.float64) * 100.0)
    return r["es"] * 0.01, r["qs"]


def entropy(tk, p, qtot, zc):
    """Parcel entropy per unit dry air (Raymond & Blyth 1992, Eq. 1).
    tk [K], p [hPa], qtot [kg/kg]."""
    L = zc["latvap"] - (zc["cpliq"] - zc["cpwv"]) * (tk - zc["tfreez"])
    _est, qst = qsat_hpa(tk, p)
    qv = jnp.minimum(qtot, qst)
    e = qv * p / (zc["epsilo"] + qv)
    return ((zc["cpair"] + qtot * zc["cpliq"])
            * jnp.log(tk / zc["tfreez"])
            - zc["rdair"] * jnp.log((p - e) / PREF)
            + L * qv / tk - qv * zc["rh2o"] * jnp.log(qv / qst))


def ientropy(s, p, qt, tfg, zc):
    """Invert entropy(T, p, qt) = s for T by Brent's method, exactly
    as zm_conv_util ientropy (bracket tfg +/- 10 K). Vectorized with
    frozen-on-convergence lanes. Returns (T, qst)."""
    s = jnp.asarray(s, dtype=jnp.float64)
    p = jnp.asarray(p, dtype=jnp.float64)
    qt = jnp.asarray(qt, dtype=jnp.float64)
    tfg = jnp.asarray(tfg, dtype=jnp.float64)

    a = tfg - 10.0
    b = tfg + 10.0
    fa = entropy(a, p, qt, zc) - s
    fb = entropy(b, p, qt, zc) - s
    cc = b
    fc = fb
    d = jnp.zeros_like(b)
    ebr = jnp.zeros_like(b)
    done = jnp.zeros(b.shape, dtype=bool)

    for _ in range(LOOPMAX + 1):
        act = ~done
        # re-bracket: move c to a when fb, fc share sign
        m1 = act & (((fb > 0.0) & (fc > 0.0))
                    | ((fb < 0.0) & (fc < 0.0)))
        cc = jnp.where(m1, a, cc)
        d = jnp.where(m1, b - a, d)
        ebr = jnp.where(m1, b - a, ebr)
        fc = jnp.where(m1, fa, fc)
        # keep the best estimate in b (a=b; b=c; c=a sequential swap)
        m2 = act & (jnp.abs(fc) < jnp.abs(fb))
        a, b, cc = (jnp.where(m2, b, a), jnp.where(m2, cc, b),
                    jnp.where(m2, b, cc))
        fa, fb, fc = (jnp.where(m2, fb, fa), jnp.where(m2, fc, fb),
                      jnp.where(m2, fb, fc))

        tol = 2.0 * TOL_EPS * jnp.abs(b) + 0.5 * TOL_COEFF
        xm = 0.5 * (cc - b)
        conv = act & ((jnp.abs(xm) <= tol) | (fb == 0.0))
        done = done | conv
        act = act & ~conv
        if not bool(jnp.any(act)):
            break

        # inverse quadratic / secant candidate step
        use_i = (jnp.abs(ebr) >= tol) & (jnp.abs(fa) > jnp.abs(fb))
        sbr = fb / fa
        secant = a == cc
        pbr = jnp.where(secant, 2.0 * xm * sbr,
                        (fb / fa) * (2.0 * xm * (fa / fc)
                                     * ((fa / fc) - (fb / fc))
                                     - (b - a) * ((fb / fc) - 1.0)))
        qbr = jnp.where(secant, 1.0 - sbr,
                        ((fa / fc) - 1.0) * ((fb / fc) - 1.0)
                        * (sbr - 1.0))
        qbr = jnp.where(pbr > 0.0, -qbr, qbr)
        pbr = jnp.abs(pbr)
        accept = 2.0 * pbr < jnp.minimum(
            3.0 * xm * qbr - jnp.abs(tol * qbr), jnp.abs(ebr * qbr))
        take = use_i & accept
        d_new = jnp.where(take, pbr / qbr, xm)
        ebr_new = jnp.where(take, d, xm)
        d = jnp.where(act, d_new, d)
        ebr = jnp.where(act, ebr_new, ebr)
        # b <- b + max(|d|, tolerance)-limited step; refresh fb
        a = jnp.where(act, b, a)
        fa = jnp.where(act, fb, fa)
        bstep = jnp.where(jnp.abs(d) > tol, d, _sign(tol, xm))
        b = jnp.where(act, b + bstep, b)
        fb = jnp.where(act, entropy(b, p, qt, zc) - s, fb)

    if bool(jnp.any(~done)):
        raise RuntimeError(
            "ientropy: Brent iteration did not converge "
            "(Fortran would endrun here)")
    t = b
    _est, qst = qsat_hpa(t, p)
    return t, qst


# ---------------------------------------------------------------------------
# find_mse_max
# ---------------------------------------------------------------------------
def find_mse_max(temperature, zmid, sp_humidity, msemax_top_k,
                 msemax_klev, num_msg, zc, zp, pergro_active=False):
    """Level of max moist static energy for parcel launch. All level
    indices 0-based. Returns (msemax_klev, mse_max_val)."""
    ncol, pver = temperature.shape
    mse_max_val = jnp.zeros(ncol)
    bot = pver - zp["mx_bot_lyr_adj"] - 1      # 0-based bottom of search
    for j in range(bot, num_msg - 1, -1):
        mse_env = (zc["cpair"] * temperature[:, j]
                   + zc["grav"] * zmid[:, j]
                   + zc["latvap"] * sp_humidity[:, j])
        if pergro_active:
            rhd = (mse_env - mse_max_val) / (mse_env + mse_max_val)
            upd = (j >= msemax_top_k) & (rhd > PERGRO_RHD_THRESHOLD)
        else:
            upd = (j >= msemax_top_k) & (mse_env > mse_max_val)
        mse_max_val = jnp.where(upd, mse_env, mse_max_val)
        msemax_klev = jnp.where(upd, j, msemax_klev)
    return msemax_klev, mse_max_val


# ---------------------------------------------------------------------------
# compute_dilute_parcel
# ---------------------------------------------------------------------------
def compute_dilute_parcel(pmid, temperature, sp_humidity, tpert, pblt,
                          klaunch, num_msg, zc, zp, parcel_temp,
                          parcel_vtemp, parcel_qsat, lcl_pmid,
                          lcl_temperature, lcl_klev):
    """Entraining parcel ascent from klaunch (0-based). Returns
    (parcel_temp, parcel_vtemp, parcel_qsat, lcl_pmid,
    lcl_temperature, lcl_klev)."""
    ncol, pver = temperature.shape

    # tpert only perturbs PBL-rooted parcels when tpert_fix is on
    if zp["tpert_fix"]:
        tpert_loc = jnp.where(klaunch < pblt, 0.0, tpert)
    else:
        tpert_loc = jnp.asarray(tpert, dtype=jnp.float64)

    tmix = jnp.zeros((ncol, pver))
    qtmix = jnp.zeros((ncol, pver))
    qsmix = jnp.zeros((ncol, pver))
    smix = jnp.zeros((ncol, pver))
    sp0 = jnp.zeros(ncol)
    qtp0 = jnp.zeros(ncol)
    mp0 = jnp.zeros(ncol)
    sp_a = jnp.zeros(ncol)
    qtp_a = jnp.zeros(ncol)
    mp_a = jnp.zeros(ncol)

    # ---- entrainment loop (sequential, surface -> top) ----
    for j in range(pver - 1, num_msg - 1, -1):
        launch = klaunch == j
        ascend = j < klaunch
        active = launch | ascend
        if not bool(jnp.any(active)):
            continue
        t_j = temperature[:, j]
        p_j = pmid[:, j]
        q_j = sp_humidity[:, j]
        # environment entropy at (t_j, p_j, q_j): this IS sp0 for the
        # launch lanes, and a safe self-consistent dummy for the rest
        s_env_j = entropy(t_j, p_j, q_j, zc)

        mp0 = jnp.where(launch, 1.0, mp0)
        qtp0 = jnp.where(launch, q_j, qtp0)
        sp0 = jnp.where(launch, s_env_j, sp0)

        any_ascend = bool(jnp.any(ascend))
        if any_ascend:
            dp = pmid[:, j] - pmid[:, j + 1]
            qtenv = 0.5 * (sp_humidity[:, j] + sp_humidity[:, j + 1])
            tenv = 0.5 * (temperature[:, j] + temperature[:, j + 1])
            penv = 0.5 * (pmid[:, j] + pmid[:, j + 1])
            senv = entropy(tenv, penv, qtenv, zc)
            dpdz = -(penv * zc["grav"]) / (zc["rdair"] * tenv)
            dmpdp = zp["dmpdz"] / dpdz     # dmpdz * dzdp
            sp_a = jnp.where(ascend, sp_a - dmpdp * dp * senv, sp_a)
            qtp_a = jnp.where(ascend, qtp_a - dmpdp * dp * qtenv,
                              qtp_a)
            mp_a = jnp.where(ascend, mp_a - dmpdp * dp, mp_a)

        smix_j = jnp.where(
            launch, sp0,
            jnp.where(ascend, (sp0 + sp_a) / (mp0 + mp_a), smix[:, j]))
        qtmix_j = jnp.where(
            launch, qtp0,
            jnp.where(ascend, (qtp0 + qtp_a) / (mp0 + mp_a),
                      qtmix[:, j]))
        smix = smix.at[:, j].set(smix_j)
        qtmix = qtmix.at[:, j].set(qtmix_j)

        jp1 = min(j + 1, pver - 1)         # tmix[:, j+1] only read
        tfg = jnp.where(launch, t_j,       # under ascend (j+1 valid)
                        jnp.where(ascend, tmix[:, jp1], t_j))
        tmix_j, qsmix_j = ientropy(
            jnp.where(active, smix_j, s_env_j), p_j,
            jnp.where(active, qtmix_j, q_j),
            jnp.where(active, tfg, t_j), zc)
        tmix = tmix.at[:, j].set(jnp.where(active, tmix_j, tmix[:, j]))
        qsmix = qsmix.at[:, j].set(
            jnp.where(active, qsmix_j, qsmix[:, j]))

        # LCL: first level with qsmix <= qtmix on ascent
        if any_ascend:
            lclm = (ascend & (qsmix[:, j] <= qtmix[:, j])
                    & (qsmix[:, j + 1] > qtmix[:, j + 1]))
            if bool(jnp.any(lclm)):
                qxsk = qtmix[:, j] - qsmix[:, j]
                qxskp1 = qtmix[:, j + 1] - qsmix[:, j + 1]
                dqxsdp = (qxsk - qxskp1) / dp
                lcl_p_new = pmid[:, j + 1] - qxskp1 / dqxsdp
                dsdp = (smix[:, j] - smix[:, j + 1]) / dp
                dqtdp = (qtmix[:, j] - qtmix[:, j + 1]) / dp
                slcl = (smix[:, j + 1]
                        + dsdp * (lcl_p_new - pmid[:, j + 1]))
                qtlcl = (qtmix[:, j + 1]
                         + dqtdp * (lcl_p_new - pmid[:, j + 1]))
                t_lcl, _qs_lcl = ientropy(
                    jnp.where(lclm, slcl, s_env_j),
                    jnp.where(lclm, lcl_p_new, p_j),
                    jnp.where(lclm, qtlcl, q_j),
                    jnp.where(lclm, tmix[:, j], t_j), zc)
                lcl_klev = jnp.where(lclm, j, lcl_klev)
                lcl_pmid = jnp.where(lclm, lcl_p_new, lcl_pmid)
                lcl_temperature = jnp.where(lclm, t_lcl,
                                            lcl_temperature)

    # ---- condensation / freezing adjustment loop ----
    xsh2o = jnp.zeros((ncol, pver))
    ds_xsh2o = jnp.zeros((ncol, pver))
    ds_freeze = jnp.zeros((ncol, pver))
    for j in range(pver - 1, num_msg - 1, -1):
        launch = klaunch == j
        ascend = j < klaunch
        if bool(jnp.any(launch)):
            pq_l = sp_humidity[:, j]
            pv_l = ((tmix[:, j] + zp["tpert_fac"] * tpert_loc)
                    * (1.0 + zc["zvir"] * pq_l) / (1.0 + pq_l))
            parcel_temp = parcel_temp.at[:, j].set(
                jnp.where(launch, tmix[:, j], parcel_temp[:, j]))
            parcel_qsat = parcel_qsat.at[:, j].set(
                jnp.where(launch, pq_l, parcel_qsat[:, j]))
            parcel_vtemp = parcel_vtemp.at[:, j].set(
                jnp.where(launch, pv_l, parcel_vtemp[:, j]))
        if not bool(jnp.any(ascend)):
            continue
        t_j = temperature[:, j]
        p_j = pmid[:, j]
        q_j = sp_humidity[:, j]
        s_env_j = entropy(t_j, p_j, q_j, zc)   # safe dummy inputs
        new_q = jnp.zeros(ncol)
        for _ii in range(NIT_LHEAT):
            xs_j = jnp.maximum(0.0,
                               qtmix[:, j] - qsmix[:, j] - LWMAX)
            dsx_j = (ds_xsh2o[:, j + 1]
                     - zc["cpliq"]
                     * jnp.log(tmix[:, j] / zc["tfreez"])
                     * jnp.maximum(0.0, xs_j - xsh2o[:, j + 1]))
            frz = tmix[:, j] <= zc["tfreez"]
            dsf_j = ds_freeze[:, j]
            one_off = frz & (ds_freeze[:, j + 1] == 0.0)
            dsf_j = jnp.where(
                one_off,
                (zc["latice"] / tmix[:, j])
                * jnp.maximum(0.0, qtmix[:, j] - qsmix[:, j] - xs_j),
                dsf_j)
            cont = frz & (ds_freeze[:, j + 1] != 0.0)
            dsf_j = jnp.where(
                cont,
                ds_freeze[:, j + 1] + (zc["latice"] / tmix[:, j])
                * jnp.maximum(0.0, qsmix[:, j + 1] - qsmix[:, j]),
                dsf_j)
            new_s = smix[:, j] + dsx_j + dsf_j
            new_q = qtmix[:, j] - xs_j
            tm_j, qs_j = ientropy(
                jnp.where(ascend, new_s, s_env_j), p_j,
                jnp.where(ascend, new_q, q_j),
                jnp.where(ascend, tmix[:, j], t_j), zc)
            tmix = tmix.at[:, j].set(
                jnp.where(ascend, tm_j, tmix[:, j]))
            qsmix = qsmix.at[:, j].set(
                jnp.where(ascend, qs_j, qsmix[:, j]))
            xsh2o = xsh2o.at[:, j].set(
                jnp.where(ascend, xs_j, xsh2o[:, j]))
            ds_xsh2o = ds_xsh2o.at[:, j].set(
                jnp.where(ascend, dsx_j, ds_xsh2o[:, j]))
            ds_freeze = ds_freeze.at[:, j].set(
                jnp.where(ascend, dsf_j, ds_freeze[:, j]))
        # parcel values: condensate reduces buoyancy when new_q > qsmix
        pq = jnp.where(new_q > qsmix[:, j], qsmix[:, j], new_q)
        pv = ((tmix[:, j] + zp["tpert_fac"] * tpert_loc)
              * (1.0 + zc["zvir"] * pq) / (1.0 + new_q))
        parcel_temp = parcel_temp.at[:, j].set(
            jnp.where(ascend, tmix[:, j], parcel_temp[:, j]))
        parcel_qsat = parcel_qsat.at[:, j].set(
            jnp.where(ascend, pq, parcel_qsat[:, j]))
        parcel_vtemp = parcel_vtemp.at[:, j].set(
            jnp.where(ascend, pv, parcel_vtemp[:, j]))

    return (parcel_temp, parcel_vtemp, parcel_qsat, lcl_pmid,
            lcl_temperature, lcl_klev)


# ---------------------------------------------------------------------------
# compute_cape_from_parcel
# ---------------------------------------------------------------------------
def compute_cape_from_parcel(temperature, tv, sp_humidity, pint,
                             msemax_klev, lcl_pmid, lcl_klev, num_msg,
                             zc, zp, parcel_qsat, parcel_temp,
                             parcel_vtemp):
    """Buoyancy integral -> CAPE and equilibrium level. Returns
    (parcel_qsat, parcel_temp, parcel_vtemp, eql_klev, cape)."""
    ncol, pver = temperature.shape
    num_cin = zp["num_cin"]
    lclok = lcl_pmid >= LCL_PRESSURE_THRESHOLD

    # buoyancy from launch to top-allowed level; outside the plume the
    # parcel arrays are reset to the environment (only for levels the
    # Fortran k-loop visits, i.e. j >= num_msg)
    ki = jnp.arange(pver)[None, :]
    inloop = ki >= num_msg
    inplume = (ki <= msemax_klev[:, None]) & lclok[:, None]
    buoy = jnp.where(inloop & inplume,
                     parcel_vtemp - tv + zp["tiedke_add"], 0.0)
    reset = inloop & ~inplume
    parcel_qsat = jnp.where(reset, sp_humidity, parcel_qsat)
    parcel_temp = jnp.where(reset, temperature, parcel_temp)
    parcel_vtemp = jnp.where(reset, tv, parcel_vtemp)

    # equilibrium-level candidates at sign changes of buoyancy,
    # keeping up to num_cin of them (counter saturates).
    # Fortran k runs to pver but the k < lcl_klev guard makes the last
    # midlevel a no-op, hence range(..., pver-1) here.
    eql_tmp = jnp.full((ncol, num_cin), pver - 1, dtype=jnp.int64)
    cnt = jnp.zeros(ncol, dtype=jnp.int64)
    for j in range(num_msg + 1, pver - 1):
        m = ((j < lcl_klev) & lclok
             & (buoy[:, j + 1] > 0.0) & (buoy[:, j] <= 0.0))
        cnt_new = jnp.minimum(num_cin, cnt + 1)
        for n in range(num_cin):
            eql_tmp = eql_tmp.at[:, n].set(
                jnp.where(m & (cnt_new == n + 1), j, eql_tmp[:, n]))
        cnt = jnp.where(m, cnt_new, cnt)

    # integrate buoyancy in the same ascending-k order as the Fortran
    dlnp = jnp.log(pint[:, 1:] / pint[:, :-1])
    cape_tmp = jnp.zeros((ncol, num_cin))
    for n in range(num_cin):
        acc = jnp.zeros(ncol)
        for j in range(num_msg, pver):
            m = (lclok & (j <= msemax_klev) & (j > eql_tmp[:, n]))
            acc = acc + jnp.where(
                m, zc["rdair"] * buoy[:, j] * dlnp[:, j], 0.0)
        cape_tmp = cape_tmp.at[:, n].set(acc)

    # keep the largest candidate CAPE (strict >, first wins ties)
    cape = jnp.zeros(ncol)
    eql_klev = jnp.full(ncol, pver - 1, dtype=jnp.int64)
    for n in range(num_cin):
        better = cape_tmp[:, n] > cape
        cape = jnp.where(better, cape_tmp[:, n], cape)
        eql_klev = jnp.where(better, eql_tmp[:, n], eql_klev)
    cape = jnp.maximum(cape, 0.0)
    return parcel_qsat, parcel_temp, parcel_vtemp, eql_klev, cape


# ---------------------------------------------------------------------------
# compute_dilute_cape (top level)
# ---------------------------------------------------------------------------
def compute_dilute_cape(sp_humidity, temperature, zmid, pmid, pint,
                        pblt, tpert, num_msg, zc, zp,
                        calc_msemax_klev=True, prev_msemax_klev=None,
                        use_input_tq_mx=False, q_mx=None, t_mx=None,
                        pergro_active=False):
    """Dilute-parcel CAPE. Inputs are (ncol, pver) C-ordered arrays,
    level 0 = top; pmid/pint in hPa; pblt/prev_msemax_klev 0-based.
    Returns a dict with parcel_temp, parcel_qsat, msemax_klev,
    lcl_temperature, lcl_klev, eql_klev, cape, q_mx, t_mx
    (level indices 0-based)."""
    q = jnp.asarray(sp_humidity, dtype=jnp.float64)
    t = jnp.asarray(temperature, dtype=jnp.float64)
    zmid = jnp.asarray(zmid, dtype=jnp.float64)
    pmid = jnp.asarray(pmid, dtype=jnp.float64)
    pint = jnp.asarray(pint, dtype=jnp.float64)
    pblt = jnp.asarray(pblt, dtype=jnp.int64)
    tpert = jnp.asarray(tpert, dtype=jnp.float64)
    ncol, pver = t.shape
    rows = jnp.arange(ncol)

    if use_input_tq_mx:
        if prev_msemax_klev is None or q_mx is None or t_mx is None:
            raise ValueError("use_input_tq_mx requires "
                             "prev_msemax_klev, q_mx and t_mx")
        msemax_klev = jnp.asarray(prev_msemax_klev, dtype=jnp.int64)
        q = q.at[rows, msemax_klev].set(jnp.asarray(q_mx))
        t = t.at[rows, msemax_klev].set(jnp.asarray(t_mx))
    else:
        msemax_klev = jnp.full(ncol, pver - 1, dtype=jnp.int64)

    tv = t * (1.0 + zc["zvir"] * q) / (1.0 + q)
    parcel_temp = t
    parcel_qsat = q
    parcel_vtemp = tv

    # unrestricted launch level (ULL): highest level below 600 hPa
    if zp["trig_ull"]:
        pblt_ull = jnp.zeros(ncol, dtype=jnp.int64)   # Fortran init 1
        for j in range(pver - 2, num_msg - 1, -1):
            m = ((pmid[:, j] <= ULL_UPPER_LAUNCH_PRESSURE)
                 & (pmid[:, j + 1] > ULL_UPPER_LAUNCH_PRESSURE))
            pblt_ull = jnp.where(m, j, pblt_ull)

    # launch level = max-MSE level (or the previous call's)
    if zp["trig_dcape"] and not calc_msemax_klev:
        if prev_msemax_klev is None:
            raise ValueError("trig_dcape with calc_msemax_klev=False "
                             "requires prev_msemax_klev")
        msemax_klev = jnp.asarray(prev_msemax_klev, dtype=jnp.int64)
    elif not use_input_tq_mx:
        top = pblt_ull if zp["trig_ull"] else pblt
        msemax_klev, _ = find_mse_max(t, zmid, q, top, msemax_klev,
                                      num_msg, zc, zp, pergro_active)

    # save launch-level T/q ("q_mx/t_mx always present" semantics)
    if use_input_tq_mx:
        q_mx_out = jnp.asarray(q_mx, dtype=jnp.float64)
        t_mx_out = jnp.asarray(t_mx, dtype=jnp.float64)
    else:
        q_mx_out = q[rows, msemax_klev]
        t_mx_out = t[rows, msemax_klev]

    lcl_klev = msemax_klev
    lcl_pmid = pmid[rows, msemax_klev]
    lcl_temperature = t[rows, msemax_klev]

    (parcel_temp, parcel_vtemp, parcel_qsat, lcl_pmid,
     lcl_temperature, lcl_klev) = compute_dilute_parcel(
        pmid, t, q, tpert, pblt, msemax_klev, num_msg, zc, zp,
        parcel_temp, parcel_vtemp, parcel_qsat, lcl_pmid,
        lcl_temperature, lcl_klev)

    (parcel_qsat, parcel_temp, parcel_vtemp, eql_klev,
     cape) = compute_cape_from_parcel(
        t, tv, q, pint, msemax_klev, lcl_pmid, lcl_klev, num_msg, zc,
        zp, parcel_qsat, parcel_temp, parcel_vtemp)

    return dict(parcel_temp=parcel_temp, parcel_qsat=parcel_qsat,
                parcel_vtemp=parcel_vtemp, msemax_klev=msemax_klev,
                lcl_temperature=lcl_temperature, lcl_klev=lcl_klev,
                lcl_pmid=lcl_pmid, eql_klev=eql_klev, cape=cape,
                q_mx=q_mx_out, t_mx=t_mx_out)
