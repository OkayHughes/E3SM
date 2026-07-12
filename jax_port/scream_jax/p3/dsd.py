"""Cloud and rain drop-size-distribution parameters.

Source: components/eamxx/src/physics/p3/impl/p3_dsd2_impl.hpp
(get_cloud_dsd2, get_rain_dsd2). iparam = 3 (Khairoutdinov-Kogan), so the
Seifert-Beheng dnu interpolation branch is compiled out in the C++ and
omitted here (nu stays 0).
"""

import jax.numpy as jnp
from jax.scipy.special import gammaln

from ..foundation import constants as c


def _tgamma(x):
    return jnp.exp(gammaln(x))


def get_cloud_dsd2(qc, nc, rho, context):
    """Cloud DSD parameters (Functions::get_cloud_dsd2).

    Returns (nc, mu_c, nu, lamc, cdist, cdist1); nc is updated (floored at
    NSMALL and recomputed where the lambda limiters bind).
    """
    qc = jnp.asarray(qc)
    nc = jnp.asarray(nc)
    rho = jnp.asarray(rho)

    gt = (qc >= c.QSMALL) & context

    nc_lim = jnp.where(gt, jnp.maximum(nc, c.NSMALL), 1.0)
    mu_l = 0.0005714 * (nc_lim * 1.0e-6 * rho) + 0.2714
    mu_l = 1.0 / (mu_l * mu_l) - 1.0
    mu_c = jnp.where(gt, jnp.clip(mu_l, 2.0, 15.0), 0.0)

    qc_safe = jnp.where(gt, qc, 1.0)
    lamc = jnp.where(
        gt,
        jnp.cbrt(c.CONS1 * nc_lim * (mu_c + 3) * (mu_c + 2) * (mu_c + 1) / qc_safe),
        0.0)

    lammin = (mu_c + 1) * 2.5e4
    lammax = (mu_c + 1) * 1.0e6
    lt_min = (lamc < lammin) & gt
    gt_max = (lamc > lammax) & gt
    lamc = jnp.where(lt_min, lammin, jnp.where(gt_max, lammax, lamc))

    either = lt_min | gt_max
    nc_new = jnp.where(
        either,
        6.0 * lamc ** 3 * qc_safe
        / (c.Pi * c.RHO_H2O * (mu_c + 3) * (mu_c + 2) * (mu_c + 1)),
        jnp.where(gt, nc_lim, nc))

    lamc_safe = jnp.where(gt, lamc, 1.0)
    cdist = jnp.where(gt, nc_new * (mu_c + 1) / lamc_safe, 0.0)
    cdist1 = jnp.where(gt, nc_new / _tgamma(mu_c + 1), 0.0)
    nu = jnp.zeros_like(qc)  # iparam != 1
    return nc_new, mu_c, nu, lamc, cdist, cdist1


def get_rain_dsd2(qr, nr, constant_mu_rain, context):
    """Rain DSD parameters (Functions::get_rain_dsd2).

    Returns (nr, mu_r, lamr); nr is updated (floored at NSMALL where rain
    present, recomputed where the lambda limiters bind).
    """
    qr = jnp.asarray(qr)
    nr = jnp.asarray(nr)

    gt = (qr >= c.QSMALL) & context

    nr_lim = jnp.maximum(nr, c.NSMALL)
    mu_r = jnp.where(gt, constant_mu_rain, jnp.zeros_like(qr))

    mass_to_d3 = c.CONS1 * (mu_r + 3) * (mu_r + 2) * (mu_r + 1)
    qr_safe = jnp.where(gt, qr, 1.0)
    lamr = jnp.where(gt, jnp.cbrt(mass_to_d3 * nr_lim / qr_safe), 0.0)

    lammax = (mu_r + 1.0) * 1.0e5
    lammin = (mu_r + 1.0) * 500.0
    lt = gt & (lamr < lammin)
    gtm = gt & (lamr > lammax)
    lamr = jnp.where(lt, lammin, jnp.where(gtm, lammax, lamr))

    nr_new = jnp.where(gt, nr_lim, nr)
    nr_new = jnp.where(lt | gtm, lamr ** 3 * qr_safe / mass_to_d3, nr_new)
    return nr_new, mu_r, lamr
