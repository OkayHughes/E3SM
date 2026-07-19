"""Cloud fraction kernel (CldFractionFunctions::main).

Physics: binary ice cloud fraction by thresholding qi (two thresholds — the
model-internal one and a looser "for output/analysis" one), total cloud
fraction as max(ice, liquid). Source:
components/eamxx/src/physics/cld_fraction/cld_fraction_main_impl.hpp
(calc_icefrac, calc_totalfrac).

Differentiability: with the default ``smooth_width=0.0`` this is the
exact (bitwise-original) binary scheme, whose qi-gradient is
identically zero — the ice fraction is piecewise constant. Passing
``smooth_width > 0`` switches on the approximation-by-identity
surrogates (scream_jax.foundation.smoothing): the threshold becomes a
sigmoid of relative exceedance ``log(qi/thresh)``-like form
(``(qi - thresh)/thresh`` here) and the max becomes smooth_max, giving
finite nonzero gradients. The surrogate converges to the binary scheme
as smooth_width -> 0.
"""

from functools import partial

import jax
import jax.numpy as jnp

from ..foundation import smoothing


@partial(jax.jit, static_argnames=("smooth_width",))
def cld_fraction_main(ice_threshold, ice_4out_threshold, qi, liq_cld_frac,
                      smooth_width=0.0):
    """Compute cloud fractions from ice mixing ratio and liquid cloud fraction.

    Args:
        ice_threshold: qi threshold [kg/kg] for the model-internal ice cloud
            fraction (strict > comparison, as in the C++).
        ice_4out_threshold: looser threshold for the `_4out` analysis fields.
        qi: ice mass mixing ratio, any shape (typically (ncol, nlev)).
        liq_cld_frac: liquid cloud fraction, same shape.
        smooth_width: 0.0 (default) = exact binary scheme; > 0 = smooth
            surrogate with transition width ``smooth_width * threshold``
            in qi (static under jit).

    Returns:
        (ice_cld_frac, tot_cld_frac, ice_cld_frac_4out, tot_cld_frac_4out)
    """
    qi = jnp.asarray(qi)
    liq_cld_frac = jnp.asarray(liq_cld_frac)

    one = jnp.ones((), dtype=liq_cld_frac.dtype)
    zero = jnp.zeros((), dtype=liq_cld_frac.dtype)

    if smooth_width == 0.0:
        ice_cld_frac = jnp.where(qi > ice_threshold, one, zero)
        ice_cld_frac_4out = jnp.where(qi > ice_4out_threshold, one, zero)
        tot_cld_frac = jnp.maximum(ice_cld_frac, liq_cld_frac)
        tot_cld_frac_4out = jnp.maximum(ice_cld_frac_4out, liq_cld_frac)
    else:
        ice_cld_frac = smoothing.step(qi - ice_threshold, smooth_width,
                                      scale=ice_threshold)
        ice_cld_frac_4out = smoothing.step(qi - ice_4out_threshold,
                                           smooth_width,
                                           scale=ice_4out_threshold)
        tot_cld_frac = smoothing.smooth_max(ice_cld_frac, liq_cld_frac,
                                            smooth_width)
        tot_cld_frac_4out = smoothing.smooth_max(ice_cld_frac_4out,
                                                 liq_cld_frac, smooth_width)

    return ice_cld_frac, tot_cld_frac, ice_cld_frac_4out, tot_cld_frac_4out
