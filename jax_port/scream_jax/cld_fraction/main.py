"""Cloud fraction kernel (CldFractionFunctions::main).

Physics: binary ice cloud fraction by thresholding qi (two thresholds — the
model-internal one and a looser "for output/analysis" one), total cloud
fraction as max(ice, liquid). Source:
components/eamxx/src/physics/cld_fraction/cld_fraction_main_impl.hpp
(calc_icefrac, calc_totalfrac).
"""

import jax
import jax.numpy as jnp


@jax.jit
def cld_fraction_main(ice_threshold, ice_4out_threshold, qi, liq_cld_frac):
    """Compute cloud fractions from ice mixing ratio and liquid cloud fraction.

    Args:
        ice_threshold: qi threshold [kg/kg] for the model-internal ice cloud
            fraction (strict > comparison, as in the C++).
        ice_4out_threshold: looser threshold for the `_4out` analysis fields.
        qi: ice mass mixing ratio, any shape (typically (ncol, nlev)).
        liq_cld_frac: liquid cloud fraction, same shape.

    Returns:
        (ice_cld_frac, tot_cld_frac, ice_cld_frac_4out, tot_cld_frac_4out)
    """
    qi = jnp.asarray(qi)
    liq_cld_frac = jnp.asarray(liq_cld_frac)

    one = jnp.ones((), dtype=liq_cld_frac.dtype)
    zero = jnp.zeros((), dtype=liq_cld_frac.dtype)

    ice_cld_frac = jnp.where(qi > ice_threshold, one, zero)
    ice_cld_frac_4out = jnp.where(qi > ice_4out_threshold, one, zero)

    tot_cld_frac = jnp.maximum(ice_cld_frac, liq_cld_frac)
    tot_cld_frac_4out = jnp.maximum(ice_cld_frac_4out, liq_cld_frac)

    return ice_cld_frac, tot_cld_frac, ice_cld_frac_4out, tot_cld_frac_4out
