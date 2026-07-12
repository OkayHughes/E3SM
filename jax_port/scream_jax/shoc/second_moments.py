"""Second-moment building blocks (fluxes, variances/covariances).

Sources (components/eamxx/src/physics/shoc/impl/):
  shoc_calc_shoc_vertflux_impl.hpp, shoc_calc_shoc_varorcovar_impl.hpp
  (the diag_second_* driver chain follows in this module.)

These kernels write only the interior interface levels k = 1..nlev-1;
boundary values are owned by the lbycond/ubycond kernels. They therefore
take the current output array and return a functional update of its
interior band.
"""

import jax
import jax.numpy as jnp


@jax.jit
def calc_shoc_vertflux(tkh_zi, dz_zi, invar, vertflux):
    """Downgradient vertical turbulent flux on interior interfaces
    (Functions::calc_shoc_vertflux):

        vertflux(k) = -tkh_zi(k) * (invar(k-1) - invar(k)) / dz_zi(k)

    invar has nlev entries; tkh_zi, dz_zi, vertflux have nlev+1. Only
    k = 1..nlev-1 of vertflux is replaced.
    """
    tkh_zi = jnp.asarray(tkh_zi)
    dz_zi = jnp.asarray(dz_zi)
    invar = jnp.asarray(invar)
    vertflux = jnp.asarray(vertflux)

    band = slice(1, invar.shape[-1])  # 1..nlev-1 inclusive
    flux = -(tkh_zi[..., band] / dz_zi[..., band]
             * (invar[..., :-1] - invar[..., 1:]))
    return vertflux.at[..., band].set(flux)


@jax.jit
def calc_shoc_varorcovar(tunefac, isotropy_zi, tkh_zi, dz_zi,
                         invar1, invar2, varorcovar):
    """SGS variance/covariance on interior interfaces
    (Functions::calc_shoc_varorcovar):

        varorcovar(k) = tunefac * isotropy_zi(k) * tkh_zi(k)
                        * (d invar1)(d invar2) / dz_zi(k)^2

    invar1/2 have nlev entries; the interface arrays nlev+1. Only
    k = 1..nlev-1 of varorcovar is replaced.
    """
    isotropy_zi = jnp.asarray(isotropy_zi)
    tkh_zi = jnp.asarray(tkh_zi)
    dz_zi = jnp.asarray(dz_zi)
    invar1 = jnp.asarray(invar1)
    invar2 = jnp.asarray(invar2)
    varorcovar = jnp.asarray(varorcovar)

    band = slice(1, invar1.shape[-1])
    d1 = invar1[..., :-1] - invar1[..., 1:]
    d2 = invar2[..., :-1] - invar2[..., 1:]
    val = (tunefac * (isotropy_zi[..., band] * tkh_zi[..., band])
           * d1 * d2 / dz_zi[..., band] ** 2)
    return varorcovar.at[..., band].set(val)
