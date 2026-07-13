"""CLUBB parcel-based mixing length — port of mixing_length.F90
(compute_mixing_length) for the EAM configuration.

PORT_NOTES
----------
Source: components/eam/src/physics/clubb/mixing_length.F90 (the
Huebler-optimized single-pass algorithm, clubb:ticket:834).

Physics: for every thermodynamic level i, an entraining parcel is
given TKE = zm2zt(em)(i) and launched upward (and, separately,
downward).  Between adjacent levels its thl/rt evolve by the exact
solution of d(x_par)/dz = -mu * (x_par - x_env) with x_env linear in
z, its buoyancy dCAPE/dz = g*(thv_par - thvm)/thvm is accumulated by
the trapezoidal rule, and the parcel stops where the (signed) CAPE
increments have exhausted its initial TKE; the sub-grid remainder is
solved with the quadratic formula (special-cased when the two
endpoint dCAPE/dz values are equal to within relative eps=1e-10).
Lscale_up/down are made nonlocal (a lower parcel that rises above a
higher one drags the higher one's Lscale_up along — the sequential
max-altitude smoothing), floored by a surface-layer minimum lminh
(lmin at ground, tapering to 0 at 500 m above gr%zm(1) under
l_implemented=T), and combined as Lscale = sqrt(up*down), with
boundary copies and the min(Lscale, Lscale_max) host-grid cap.

Faithfulness details (all verbatim):
  - The initial one-level ascent/descent is precomputed VECTORIZED
    over all levels (thl_par_1/rt_par_1/dCAPE_dz_1/CAPE_incr_1) and
    the per-parcel loop resumes from level i+-2 — the port reproduces
    this exact structure and FP order, including
    Lscale_up(i) = ((zlmin + zt(j-1)) - zt(i)) association.
  - compute_rsat_parcel: l_include_ice=.false. (compile-time) and
    l_sat_mixrat_lookup=.false. (model_flags default) reduce it to
    sat_mixrat_liq (flatau); rsatl_par = liq + (ice-liq)*0 evaluates
    to exactly liq.
  - The mu==0 abort guard is not ported (EAMv3 mu=5e-4; callers pass
    a positive entrainment rate).
  - length_check / stats are debug/stats-only (debug_level 0,
    l_stats=F): not ported.
  - Level 1 (the below-ground ghost): Lscale_up = Lscale_down = 0,
    Lscale(1) = Lscale(2); Lscale(nz) = Lscale(nz-1).
  - Lscale_up/Lscale_down keep zlmin = 0.1 m at the levels the main
    loops do not visit (up: nz-1, nz; down: 2), subject to the lminh
    floor.

DEAD in EAM (documented, not ported): the l_avg_Lscale perturbed-
Lscale averaging in advance_clubb_core is gated by a COMPILE-TIME
`logical, parameter :: l_avg_Lscale = .false.` (and
l_Lscale_plume_centered = .false. in model_flags), so
compute_mixing_length is called exactly ONCE per timestep with the
unperturbed thlm/rtm and mu = the namelist tunable (newmu = mu;
CLUBBND_CAM is not defined).  Lscale_pert_1/2 are unused_var.

The mu argument must be passed as a traced scalar (the segment driver
wraps it in jnp.asarray) so XLA cannot constant-fold divisions by mu
into reciprocal multiplies.
"""

import jax.numpy as jnp
from jax import lax

from .constants import (CP, EP, EP1, EP2, EPS, GRAV, LV, RD,
                        ZERO_THRESHOLD)
from .grid import Grid, zm2zt
from .saturation import sat_mixrat_liq

ZLMIN = 0.1                    # minimum Lscale_up/down seed [m]
LSCALE_SFCLYR_DEPTH = 500.0    # surface-layer lminh taper depth [m]


def compute_mixing_length(gr: Grid, thvm, thlm, rtm, em, lscale_max,
                          p_in_pa, exner, thv_ds, mu, lmin):
    """Port of compute_mixing_length (l_implemented=.true.).

    All profile arguments are (nz,) on their Fortran grids (thvm/thlm/
    rtm/exner/p/thv_ds on zt, em on zm); lscale_max, mu, lmin scalars.
    Returns (Lscale, Lscale_up, Lscale_down).
    """
    thvm = jnp.asarray(thvm)
    nz = thvm.shape[0]
    mu = jnp.asarray(mu, dtype=jnp.float64)

    # ---- precalculations (verbatim order) ---------------------------
    exp_mu_dzm = jnp.exp(-mu * gr.dzm)
    invrs_dzm_on_mu = gr.invrs_dzm / mu
    grav_on_thvm = GRAV / thvm
    lv_coef = LV / (exner * CP) - EP2 * thv_ds
    entrain_coef = (1.0 - exp_mu_dzm) * invrs_dzm_on_mu
    lv2_coef = EP * LV ** 2 / (RD * CP)
    invrs_sfclyr_depth = 1.0 / LSCALE_SFCLYR_DEPTH

    tke_i = zm2zt(gr, em)
    zt = gr.zt
    dzm = gr.dzm
    invrs_dzm = gr.invrs_dzm

    def parcel_thermo(j, thl_par, rt_par):
        """s/rc/thv/dCAPE of a parcel with (thl_par, rt_par) at level
        j (the shared body of the initial and continued steps)."""
        tl_par = thl_par * exner[j]
        rsatl_par = sat_mixrat_liq(p_in_pa[j], tl_par)
        tl_sqd = tl_par ** 2
        s_par = (rt_par - rsatl_par) * tl_sqd \
            / (tl_sqd + lv2_coef * rsatl_par)
        rc_par = jnp.maximum(s_par, ZERO_THRESHOLD)
        thv_par = thl_par + EP1 * thv_ds[j] * rt_par \
            + lv_coef[j] * rc_par
        return grav_on_thvm[j] * (thv_par - thvm[j])

    # ================= upward length scale ===========================
    # recursion precalc (Fortran j = 2..nz-1; py index j uses [j-1])
    z1 = jnp.zeros(1)
    tp_thl_up = jnp.concatenate([z1, thlm[1:] - thlm[:-1]
                                 * exp_mu_dzm[:-1]
                                 - (thlm[1:] - thlm[:-1])
                                 * entrain_coef[:-1]])
    tp_rt_up = jnp.concatenate([z1, rtm[1:] - rtm[:-1]
                                * exp_mu_dzm[:-1]
                                - (rtm[1:] - rtm[:-1])
                                * entrain_coef[:-1]])

    # initial one-level ascent (Fortran j = 3..nz; vectorized j >= 1)
    thl1_up = jnp.concatenate([thlm[:1], thlm[1:] - (thlm[1:]
                               - thlm[:-1]) * entrain_coef[:-1]])
    tl1_up = thl1_up * exner
    rt1_up = jnp.concatenate([rtm[:1], rtm[1:] - (rtm[1:] - rtm[:-1])
                              * entrain_coef[:-1]])
    rsatl1_up = sat_mixrat_liq(p_in_pa, tl1_up)
    tl1_sqd = tl1_up ** 2
    s1_up = (rt1_up - rsatl1_up) * tl1_sqd \
        / (tl1_sqd + lv2_coef * rsatl1_up)
    rc1_up = jnp.maximum(s1_up, ZERO_THRESHOLD)
    thv1_up = thl1_up + EP1 * thv_ds * rt1_up + lv_coef * rc1_up
    dcape1_up = grav_on_thvm * (thv1_up - thvm)
    cape_incr1_up = jnp.concatenate(
        [z1, 0.5 * dcape1_up[1:] * dzm[:-1]])

    def up_body(i, carry):
        lup_arr, max_alt = carry
        entered = tke_i[i] + cape_incr1_up[i + 1] > 0.0

        # -- branch A: the parcel survives the first level ------------
        def while_cond(st):
            j, _thl, _rt, _dcp, _tke, exhausted, _dcj = st
            return (j < nz - 1) & jnp.logical_not(exhausted)

        def while_body(st):
            j, thl_par, rt_par, dcape_prev, tke, _exh, _dcj = st
            thl_par = tp_thl_up[j] + thl_par * exp_mu_dzm[j - 1]
            rt_par = tp_rt_up[j] + rt_par * exp_mu_dzm[j - 1]
            dcape_j = parcel_thermo(j, thl_par, rt_par)
            cape_incr = 0.5 * (dcape_j + dcape_prev) * dzm[j - 1]
            stop = tke + cape_incr <= 0.0
            return (jnp.where(stop, j, j + 1), thl_par, rt_par,
                    jnp.where(stop, dcape_prev, dcape_j),
                    jnp.where(stop, tke, tke + cape_incr),
                    stop, dcape_j)

        st0 = (i + 2, thl1_up[i + 1], rt1_up[i + 1], dcape1_up[i + 1],
               tke_i[i] + cape_incr1_up[i + 1], jnp.asarray(False),
               jnp.asarray(0.0))
        j, _thl, _rt, dcape_prev, tke, exhausted, dcape_j = \
            lax.while_loop(while_cond, while_body, st0)
        lup_a = (ZLMIN + zt[j - 1]) - zt[i]
        # partial sub-grid distance when TKE was exhausted (j < nz)
        equal = jnp.abs(dcape_j - dcape_prev) * 2 \
            <= jnp.abs(dcape_j + dcape_prev) * EPS
        part_eq = lup_a + (-tke / dcape_j)
        invrs_diff = 1.0 / (dcape_j - dcape_prev)
        part_ne = lup_a - dcape_prev * invrs_diff * dzm[j - 1] \
            - jnp.sqrt(dcape_prev ** 2
                       - 2.0 * tke * invrs_dzm[j - 1]
                       * (dcape_j - dcape_prev)) \
            * invrs_diff * dzm[j - 1]
        lup_a = jnp.where(exhausted,
                          jnp.where(equal, part_eq, part_ne), lup_a)

        # -- branch B: exhausted before one full level -----------------
        lup_b = ZLMIN - jnp.sqrt(-2.0 * tke_i[i] * dzm[i]
                                 * dcape1_up[i + 1]) / dcape1_up[i + 1]

        lup_i = jnp.where(entered, lup_a, lup_b)

        # nonlocal smoothing (sequential in i)
        below = zt[i] + lup_i < max_alt
        lup_i = jnp.where(below, max_alt - zt[i], lup_i)
        max_alt = jnp.where(below, max_alt, lup_i + zt[i])
        return lup_arr.at[i].set(lup_i), max_alt

    lscale_up = jnp.full(nz, ZLMIN).at[0].set(0.0)
    lscale_up, _ = lax.fori_loop(1, nz - 2, up_body,
                                 (lscale_up, jnp.asarray(0.0)))

    # ================= downward length scale =========================
    # recursion precalc (Fortran j = 2..nz-1; py index j uses [j])
    tp_thl_dn = jnp.concatenate([thlm[:-1] - thlm[1:]
                                 * exp_mu_dzm[:-1]
                                 - (thlm[:-1] - thlm[1:])
                                 * entrain_coef[:-1], z1])
    tp_rt_dn = jnp.concatenate([rtm[:-1] - rtm[1:] * exp_mu_dzm[:-1]
                                - (rtm[:-1] - rtm[1:])
                                * entrain_coef[:-1], z1])

    # initial one-level descent (Fortran j = 2..nz-1)
    thl1_dn = jnp.concatenate([thlm[:-1] - (thlm[:-1] - thlm[1:])
                               * entrain_coef[:-1], thlm[-1:]])
    tl1_dn = thl1_dn * exner
    rt1_dn = jnp.concatenate([rtm[:-1] - (rtm[:-1] - rtm[1:])
                              * entrain_coef[:-1], rtm[-1:]])
    rsatl1_dn = sat_mixrat_liq(p_in_pa, tl1_dn)
    tl1_dn_sqd = tl1_dn ** 2
    s1_dn = (rt1_dn - rsatl1_dn) * tl1_dn_sqd \
        / (tl1_dn_sqd + lv2_coef * rsatl1_dn)
    rc1_dn = jnp.maximum(s1_dn, ZERO_THRESHOLD)
    thv1_dn = thl1_dn + EP1 * thv_ds * rt1_dn + lv_coef * rc1_dn
    dcape1_dn = grav_on_thvm * (thv1_dn - thvm)
    cape_incr1_dn = 0.5 * dcape1_dn * dzm

    def dn_body(t, carry):
        ldn_arr, min_alt = carry
        i = nz - 1 - t
        entered = tke_i[i] - cape_incr1_dn[i - 1] > 0.0

        def while_cond(st):
            j, _thl, _rt, _dcp, _tke, exhausted, _dcj = st
            return (j >= 1) & jnp.logical_not(exhausted)

        def while_body(st):
            j, thl_par, rt_par, dcape_prev, tke, _exh, _dcj = st
            thl_par = tp_thl_dn[j] + thl_par * exp_mu_dzm[j]
            rt_par = tp_rt_dn[j] + rt_par * exp_mu_dzm[j]
            dcape_j = parcel_thermo(j, thl_par, rt_par)
            cape_incr = 0.5 * (dcape_j + dcape_prev) * dzm[j]
            stop = tke - cape_incr <= 0.0
            return (jnp.where(stop, j, j - 1), thl_par, rt_par,
                    jnp.where(stop, dcape_prev, dcape_j),
                    jnp.where(stop, tke, tke - cape_incr),
                    stop, dcape_j)

        st0 = (i - 2, thl1_dn[i - 1], rt1_dn[i - 1], dcape1_dn[i - 1],
               tke_i[i] - cape_incr1_dn[i - 1], jnp.asarray(False),
               jnp.asarray(0.0))
        j, _thl, _rt, dcape_prev, tke, exhausted, dcape_j = \
            lax.while_loop(while_cond, while_body, st0)
        ldn_a = (ZLMIN + zt[i]) - zt[j + 1]
        equal = jnp.abs(dcape_j - dcape_prev) * 2 \
            <= jnp.abs(dcape_j + dcape_prev) * EPS
        part_eq = ldn_a + (tke / dcape_j)
        invrs_diff = 1.0 / (dcape_j - dcape_prev)
        part_ne = ldn_a - dcape_prev * invrs_diff * dzm[j] \
            + jnp.sqrt(dcape_prev ** 2
                       + 2.0 * tke * invrs_dzm[j]
                       * (dcape_j - dcape_prev)) \
            * invrs_diff * dzm[j]
        ldn_a = jnp.where(exhausted,
                          jnp.where(equal, part_eq, part_ne), ldn_a)

        ldn_b = ZLMIN + jnp.sqrt(2.0 * tke_i[i] * dzm[i - 1]
                                 * dcape1_dn[i - 1]) / dcape1_dn[i - 1]

        ldn_i = jnp.where(entered, ldn_a, ldn_b)

        above = zt[i] - ldn_i > min_alt
        ldn_i = jnp.where(above, zt[i] - min_alt, ldn_i)
        min_alt = jnp.where(above, min_alt, zt[i] - ldn_i)
        return ldn_arr.at[i].set(ldn_i), min_alt

    lscale_down = jnp.full(nz, ZLMIN).at[0].set(0.0)
    lscale_down, _ = lax.fori_loop(0, nz - 2, dn_body,
                                   (lscale_down, zt[nz - 1]))

    # ================= final Lscale ==================================
    lminh = jnp.maximum(ZERO_THRESHOLD,
                        LSCALE_SFCLYR_DEPTH - (zt - gr.zm[0])) \
        * lmin * invrs_sfclyr_depth
    lscale_up = lscale_up.at[1:].set(
        jnp.maximum(lminh[1:], lscale_up[1:]))
    lscale_down = lscale_down.at[1:].set(
        jnp.maximum(lminh[1:], lscale_down[1:]))

    lscale = jnp.sqrt(lscale_up * lscale_down)
    lscale = lscale.at[0].set(lscale[1])
    lscale = lscale.at[nz - 1].set(lscale[nz - 2])
    lscale = jnp.minimum(lscale, lscale_max)

    return lscale, lscale_up, lscale_down
