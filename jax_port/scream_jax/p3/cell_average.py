"""Cell-average scaling and time/space physical variables (batch 4a).

Sources (components/eamxx/src/physics/p3/impl/):
  p3_back_to_cell_average_impl.hpp, p3_get_time_space_phys_variables_impl.hpp
"""

import functools

import jax.numpy as jnp

from ..foundation import constants as c

# Tendencies scaled by each overlap fraction in back_to_cell_average.
# (name -> which fraction), transcribed line-by-line from the C++.
_SCALE_MAP = {
    "qc2qr_accret_tend": "lr", "qr2qv_evap_tend": "r",
    "qc2qr_autoconv_tend": "l", "nc_accret_tend": "lr",
    "nc_selfcollect_tend": "l", "nc2nr_autoconv_tend": "l",
    "nr_selfcollect_tend": "r", "nr_evap_tend": "r", "ncautr": "lr",
    "qi2qv_sublim_tend": "i*", "nr_ice_shed_tend": "il",
    "qc2qi_hetero_freeze_tend": "il", "qr2qi_collect_tend": "ir",
    "qc2qr_ice_shed_tend": "il", "qi2qr_melt_tend": "i",
    "qc2qi_collect_tend": "il", "qr2qi_immers_freeze_tend": "r",
    "ni2nr_melt_tend": "i", "nc_collect_tend": "il", "ncshdc": "il",
    "nc2ni_immers_freeze_tend": "l", "nr_collect_tend": "ir",
    "ni_selfcollect_tend": "i", "qv2qi_vapdep_tend": "i*",
    "nr2ni_immers_freeze_tend": "r", "ni_sublim_tend": "i*",
    "qc2qi_berg_tend": "il",
    "ncheti_cnt": "l", "qcheti_cnt": "l", "nicnt": "l", "qicnt": "l",
    "ninuc_cnt": "l", "qinuc_cnt": "l",
}


def back_to_cell_average(cld_frac_l, cld_frac_r, cld_frac_i, tends,
                         context, use_separate_ice_liq_frac: bool = False):
    """Scale in-cloud process rates to grid-cell averages
    (Functions::back_to_cell_average).

    tends is a dict of tendency arrays keyed by the C++ names in
    _SCALE_MAP; returns a new dict with all entries scaled. Fractions:
    l/r/i = liquid/rain/ice cloud fraction, lr/il/ir = pairwise overlaps,
    'i*' = glaciated fraction when use_separate_ice_liq_frac else i.
    """
    l = jnp.asarray(cld_frac_l)
    r = jnp.asarray(cld_frac_r)
    i = jnp.asarray(cld_frac_i)
    fr = {
        "l": l, "r": r, "i": i,
        "ir": jnp.minimum(i, r), "il": jnp.minimum(i, l),
        "lr": jnp.minimum(l, r),
    }
    fr["i*"] = (jnp.maximum(0.0001, i - fr["il"])
                if use_separate_ice_liq_frac else i)

    out = {}
    for name, arr in tends.items():
        arr = jnp.asarray(arr)
        if name not in _SCALE_MAP:
            # e.g. qv2qi_nucleat_tend / ni_nucleat_tend: passed to the C++
            # function but not scaled (already cell-average rates)
            out[name] = arr
            continue
        f = fr[_SCALE_MAP[name]]
        out[name] = jnp.where(context, arr * f, arr)
    return out


def get_time_space_phys_variables(T_atm, pres, rho, qv_sat_l, qv_sat_i,
                                  context):
    """Air property and psychrometric variables
    (Functions::get_time_space_phys_variables). Returns
    (mu, dv, sc, dqsdt, dqsidt, ab, abi, kap, eii)."""
    T_atm = jnp.asarray(T_atm)
    pres = jnp.asarray(pres)

    mu = 1.496e-6 * T_atm ** 1.5 / (T_atm + 120.0)
    dv = 8.794e-5 * T_atm ** 1.81 / pres
    sc = mu / (jnp.asarray(rho) * dv)

    tval1, tval2, dtval = 253.15, 273.15, 20.0
    dum = 1.0 / (c.RV * T_atm ** 2)
    dqsdt = c.LatVap * jnp.asarray(qv_sat_l) * dum
    dqsidt = (c.LatVap + c.LatIce) * jnp.asarray(qv_sat_i) * dum
    ab = 1.0 + dqsdt * c.LatVap * c.INV_CP
    abi = 1.0 + dqsidt * (c.LatVap + c.LatIce) * c.INV_CP
    kap = 1.414e3 * mu

    eii = jnp.where(T_atm < tval1, 0.001,
                    jnp.where(T_atm < tval2,
                              0.001 + (T_atm - tval1) * (0.3 - 0.001) / dtval,
                              0.3))
    return mu, dv, sc, dqsdt, dqsidt, ab, abi, kap, eii
