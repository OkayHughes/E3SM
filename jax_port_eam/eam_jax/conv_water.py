"""Port of conv_water_4rad (eam/src/physics/cam/conv_water.F90):
grid-box average liquid and ice from stratus + cumulus, computed just
for radiation, plus the COSP shallow-condensate pbuf outputs.

PORT_NOTES
----------
- Purely per-point (i,k) arithmetic; the pbuf/phys_state plumbing of
  the Fortran becomes plain (ncol, pver) array arguments. All jnp,
  fully vectorized, no level loops.
- Switches (python-level, resolved before tracing):
  * conv_water_mode: 1 = area-weighted arithmetic average (EAMv3
    phys="default" conv_water_in_rad = 1), 2 = "arithmetic in
    emissivity" log-average. Mode 0 exists only OUTSIDE the routine
    (cloud_diagnostics.F90 skips the call), and the Fortran `case
    default` bodies are commented-out endruns that would leave
    cu_icwmr/tot_icwmr as stale loop-carried scalars (undefined
    behaviour) - the port only accepts 1 or 2.
  * zm_microp (zm_param%zm_microp, EAMv3 default zmconv_microp =
    .true.): the convective-microphysics branch partitions convective
    condensate with dp_icimr from the ZM microphysics instead of the
    stratiform ice fraction wrk1; it never reads conv_water_mode, and
    its totg_* come from the in-branch assignments because the final
    non-microp repartition block is skipped.
  * microp_scheme: 'RK' uses kabsi = 0.005 + 1/rei, anything else
    (EAMv3: 'P3') clamps rei to [13, 130] first. Only observable in
    mode 2 (kabs feeds alpha which only the emissivity average uses).
  * pergro_mods (default .false.): when the stratiform in-cloud
    condensate ls_icwmr < 100*ic_limit, wrk1 is replaced by a pure
    temperature ramp (1 below 243 K, linear to 0 at 263 K).
- Faithful oddities preserved:
  * in the "no convective cloud" branch the Fortran assigns
    tot_ice/totg_ice etc. directly, but (non-microp) the final
    repartition block OVERWRITES totg_* as tot0_frac*tot_icwmr*wrk1 -
    so e.g. a point with ast >= frac_limit gets qi replaced by
    (ql+qi)/max(frac_limit,ast)*ast*wrk1. The port reproduces the
    overwrite, not the dead in-branch stores (which survive only
    under zm_microp).
  * `sh_iclmr = sh_icwmr*(1-wrk1)` etc. keep the Fortran grouping
    sh0_frac*(sh_icwmr*wrk1) for bit agreement.
  * fice NaN guard: the FICE pbuf field can be unset in the real
    model; NaN rows produce exactly zero sh_cldliq/sh_cldice (the
    only outputs fice touches).
- Constants verbatim from conv_water.F90: kabsl = 0.090361,
  frac_limit = 0.01, ic_limit = 1e-12; gravit from shr_const_mod via
  eam_jax.constants.
- Everything float64; safe-where guards keep masked-off lanes finite
  (log/exp in mode 2, 1/rei in RK).
- Replay agreement: exact-branch arithmetic (mode 1, microp) <= 1e-13
  relative; mode 2 carries 1-ulp cross-libm exp/log noise amplified
  by the log(1+alpha*w)/alpha cancellation when the in-cloud water
  sits just above ic_limit (absolute error bounded by ulp(1)/|alpha|
  < 4e-19 kg/kg; see tests/test_conv_water.py).
"""

import jax.numpy as jnp

from .constants import GRAVIT

__all__ = ["conv_water_4rad"]

KABSL = 0.090361       # longwave liquid absorption coeff [m2/g]
FRAC_LIMIT = 0.01      # cloud-fraction significance threshold
IC_LIMIT = 1.0e-12     # in-cloud condensate significance threshold


def conv_water_4rad(t, pdel, q_cldliq, q_cldice, sh_icwmr, dp_icwmr,
                    dp_icimr, fice, sh_frac, dp_frac, ast, rei,
                    conv_water_mode=1, zm_microp=True,
                    microp_scheme="P3", pergro_mods=False):
    """Grid-box average in-cloud liquid/ice for radiation. All array
    inputs are (ncol, pver) C-ordered, level 0 = model top; t [K],
    pdel [Pa], q_* grid-box mixing ratios [kg/kg], *_icwmr/icimr
    in-cloud mixing ratios [kg/kg], fractions dimensionless, rei
    [micron], fice may contain NaN. Returns a dict with totg_liq,
    totg_ice, sh_cldliq, sh_cldice."""
    if conv_water_mode not in (1, 2):
        raise ValueError("conv_water_mode must be 1 or 2 (mode 0 "
                         "means conv_water_4rad is not called)")
    t = jnp.asarray(t, dtype=jnp.float64)
    pdel = jnp.asarray(pdel, dtype=jnp.float64)
    ql = jnp.asarray(q_cldliq, dtype=jnp.float64)
    qi = jnp.asarray(q_cldice, dtype=jnp.float64)
    sh_icwmr = jnp.asarray(sh_icwmr, dtype=jnp.float64)
    dp_icwmr = jnp.asarray(dp_icwmr, dtype=jnp.float64)
    dp_icimr = jnp.asarray(dp_icimr, dtype=jnp.float64)
    fice = jnp.asarray(fice, dtype=jnp.float64)
    sh_frac = jnp.asarray(sh_frac, dtype=jnp.float64)
    dp_frac = jnp.asarray(dp_frac, dtype=jnp.float64)
    ast = jnp.asarray(ast, dtype=jnp.float64)
    rei = jnp.asarray(rei, dtype=jnp.float64)

    sh0 = jnp.where((sh_frac <= FRAC_LIMIT) | (sh_icwmr <= IC_LIMIT),
                    0.0, sh_frac)
    dp0 = jnp.where((dp_frac <= FRAC_LIMIT) | (dp_icwmr <= IC_LIMIT),
                    0.0, dp_frac)
    cu0 = sh0 + dp0

    # stratiform condensate ice fraction (emissivity partition)
    wrk1 = jnp.minimum(1.0, jnp.maximum(0.0, qi / (qi + ql + 1.0e-36)))

    branch_a = (cu0 < FRAC_LIMIT) | ((sh_icwmr + dp_icwmr) < IC_LIMIT)

    # ---- branch A: no significant convective cloud ----
    ls_small = ast < FRAC_LIMIT
    ls_frac_a = jnp.where(ls_small, 0.0, ast)
    ls_icwmr_a = jnp.where(
        ls_small, 0.0, (ql + qi) / jnp.maximum(FRAC_LIMIT, ast))
    tot0_a = ls_frac_a
    tot_icwmr_a = ls_icwmr_a
    # in-branch grid-box stores (survive only under zm_microp)
    totg_ice_a = jnp.where(ls_small, 0.0, qi)
    totg_liq_a = jnp.where(ls_small, 0.0, ql)

    # ---- branch B: significant convective cloud ----
    cu0_s = jnp.maximum(FRAC_LIMIT, cu0)

    if zm_microp:
        conv_ice = (sh0 * (sh_icwmr * wrk1)
                    + dp0 * dp_icimr) / cu0_s
        conv_liq = (sh0 * (sh_icwmr * (1.0 - wrk1))
                    + dp0 * (dp_icwmr - dp_icimr)) / cu0_s
        tot0_b = ast + cu0
        tot0_bs = jnp.maximum(FRAC_LIMIT, tot0_b)
        tot_ice_b = (qi + cu0 * conv_ice) / tot0_bs
        tot_liq_b = (ql + cu0 * conv_liq) / tot0_bs
        totg_ice = jnp.where(branch_a, totg_ice_a, tot0_b * tot_ice_b)
        totg_liq = jnp.where(branch_a, totg_liq_a, tot0_b * tot_liq_b)
    else:
        # radiation constants for the emissivity-weighted average
        if microp_scheme == "RK":
            kabsi = 0.005 + 1.0 / jnp.where(rei != 0.0, rei, 1.0)
        else:
            kabsi = 0.005 + 1.0 / jnp.minimum(
                jnp.maximum(13.0, rei), 130.0)
        kabs = KABSL * (1.0 - wrk1) + kabsi * wrk1
        alpha = -1.66 * kabs * pdel / GRAVIT * 1000.0

        if conv_water_mode == 1:
            cu_icwmr_b = (sh0 * sh_icwmr + dp0 * dp_icwmr) / cu0_s
        else:
            sh0e = jnp.exp(alpha * sh_icwmr)
            dp0e = jnp.exp(alpha * dp_icwmr)
            arg = (sh0 * sh0e + dp0 * dp0e) / cu0_s
            cu_icwmr_b = jnp.log(jnp.where(arg > 0.0, arg, 1.0))
            cu_icwmr_b = cu_icwmr_b / alpha

        ls_icwmr_b = (ql + qi) / jnp.maximum(FRAC_LIMIT, ast)
        tot0_b = ast + cu0
        if conv_water_mode == 1:
            tot_icwmr_b = (ast * ls_icwmr_b + cu0 * cu_icwmr_b) \
                / jnp.maximum(FRAC_LIMIT, tot0_b)
        else:
            arg = (ast * jnp.exp(alpha * ls_icwmr_b)
                   + cu0 * jnp.exp(alpha * cu_icwmr_b)) \
                / jnp.maximum(FRAC_LIMIT, tot0_b)
            tot_icwmr_b = jnp.log(jnp.where(arg > 0.0, arg, 1.0))
            tot_icwmr_b = tot_icwmr_b / alpha

        ls_icwmr = jnp.where(branch_a, ls_icwmr_a, ls_icwmr_b)
        tot0 = jnp.where(branch_a, tot0_a, tot0_b)
        tot_icwmr = jnp.where(branch_a, tot_icwmr_a, tot_icwmr_b)

        # pergro: no significant stratiform condensate -> repartition
        # the convective condensate by temperature alone
        if pergro_mods:
            wrk1_t = jnp.where(
                t < 243.0, 1.0,
                jnp.where(t < 263.0, 1.0 - (t - 243.0) / 20.0, 0.0))
            wrk1 = jnp.where(ls_icwmr < 100.0 * IC_LIMIT, wrk1_t,
                             wrk1)

        # final repartition overwrites every point (both branches)
        totg_ice = tot0 * tot_icwmr * wrk1
        totg_liq = tot0 * tot_icwmr * (1.0 - wrk1)

    # COSP shallow-condensate outputs (guarded against unset FICE)
    fice_nan = jnp.isnan(fice)
    fice_s = jnp.where(fice_nan, 0.0, fice)
    sh_cldliq = jnp.where(fice_nan, 0.0,
                          sh_icwmr * (1.0 - fice_s) * sh_frac)
    sh_cldice = jnp.where(fice_nan, 0.0, sh_icwmr * fice_s * sh_frac)

    return dict(totg_liq=totg_liq, totg_ice=totg_ice,
                sh_cldliq=sh_cldliq, sh_cldice=sh_cldice)
