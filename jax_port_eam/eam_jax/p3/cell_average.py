"""Cell-average scaling and time/space physical variables.

Source: micro_p3.F90 back_to_cell_average / get_time_space_phys_variables
(eam variant). COPIED from scream_jax/p3/cell_average.py; adaptations:

  * no use_separate_ice_liq_frac branch (does not exist in EAM);
    qi2qv_sublim/qv2qi_vapdep/ni_sublim scale with cld_frac_i plainly.
  * tendency name map matches the EAM back_to_cell_average argument
    list one-to-one (including the six *_cnt terms * cld_frac_l).
"""

import jax.numpy as jnp

from . import constants as c

# tendency name -> overlap fraction, transcribed line-by-line from
# back_to_cell_average (l/r/i = liquid/rain/ice cloud fraction,
# lr/il/ir = pairwise minima). qinuc/ni_nucleat_tend are already
# cell-average and not scaled.
_SCALE_MAP = {
    "qc2qr_accret_tend": "lr", "qr2qv_evap_tend": "r",
    "qc2qr_autoconv_tend": "l", "nc_accret_tend": "lr",
    "nc_selfcollect_tend": "l", "nc2nr_autoconv_tend": "l",
    "nr_selfcollect_tend": "r", "nr_evap_tend": "r", "ncautr": "lr",
    "qi2qv_sublim_tend": "i", "nr_ice_shed_tend": "il",
    "qc2qi_hetero_freeze_tend": "il", "qrcol": "ir",
    "qc2qr_ice_shed_tend": "il", "qi2qr_melt_tend": "i",
    "qccol": "il", "qr2qi_immers_freeze_tend": "r",
    "ni2nr_melt_tend": "i", "nc_collect_tend": "il", "ncshdc": "il",
    "nc2ni_immers_freeze_tend": "l", "nr_collect_tend": "ir",
    "ni_selfcollect_tend": "i", "qidep": "i",
    "nr2ni_immers_freeze_tend": "r", "ni_sublim_tend": "i",
    "qiberg": "il",
    "ncheti_cnt": "l", "qcheti_cnt": "l", "nicnt": "l", "qicnt": "l",
    "ninuc_cnt": "l", "qinuc_cnt": "l",
}


def back_to_cell_average(cld_frac_l, cld_frac_r, cld_frac_i, tends, context):
    """Scale in-cloud process rates to grid-cell averages
    (back_to_cell_average). tends is a dict keyed per _SCALE_MAP;
    returns a new dict with all entries scaled."""
    l = jnp.asarray(cld_frac_l)
    r = jnp.asarray(cld_frac_r)
    i = jnp.asarray(cld_frac_i)
    fr = {
        "l": l, "r": r, "i": i,
        "ir": jnp.minimum(i, r), "il": jnp.minimum(i, l),
        "lr": jnp.minimum(l, r),
    }

    out = {}
    for name, arr in tends.items():
        arr = jnp.asarray(arr)
        if name not in _SCALE_MAP:
            # qinuc / ni_nucleat_tend: already cell-averaged
            out[name] = arr
            continue
        f = fr[_SCALE_MAP[name]]
        out[name] = jnp.where(context, arr * f, arr)
    return out


def get_time_space_phys_variables(T_atm, pres, rho, qv_sat_l, qv_sat_i,
                                  context):
    """Air property and psychrometric variables
    (get_time_space_phys_variables). Returns
    (mu, dv, sc, dqsdt, dqsidt, ab, abi, kap, eii)."""
    T_atm = jnp.asarray(T_atm)
    pres = jnp.asarray(pres)

    mu = 1.496e-6 * T_atm ** 1.5 / (T_atm + 120.0)
    dv = 8.794e-5 * T_atm ** 1.81 / pres
    sc = mu / (jnp.asarray(rho) * dv)

    dum = 1.0 / (c.rv * (T_atm * T_atm))
    dqsdt = c.latvap * jnp.asarray(qv_sat_l) * dum
    dqsidt = c.latsub * jnp.asarray(qv_sat_i) * dum
    ab = 1.0 + dqsdt * c.latvap * c.inv_cp
    abi = 1.0 + dqsidt * c.latsub * c.inv_cp
    kap = 1.414e3 * mu

    # simple temperature-dependent aggregation efficiency
    eii = jnp.where(T_atm < 253.15, 0.001,
                    jnp.where(T_atm < 273.15,
                              0.001 + (T_atm - 253.15) * (0.3 - 0.001) / 20.0,
                              0.3))
    return mu, dv, sc, dqsdt, dqsidt, ab, abi, kap, eii
