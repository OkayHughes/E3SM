"""Cloud and rain drop-size-distribution parameters.

Source: micro_p3.F90 get_cloud_dsd2 / get_rain_dsd2 (eam variant).
COPIED from scream_jax/p3/dsd.py; adaptations against the EAM Fortran:

  * get_rain_dsd2: mu_r = mu_r_constant (0 in EAM); lammin =
    (mu_r+1)/p3_max_mean_rain_size (namelist; scream: hard 500 = 1/2mm);
    nr recompute at the limiters via
    exp(3*log(lamr)+log(qr)+log(G(mu+1))-log(G(mu+4)))/cons1
    (scream: lamr^3*qr/mass_to_d3); cdistr = nr/G(mu_r+1) and logn0r =
    log10(nr)+(mu_r+1)*log10(lamr)-log10(G(mu_r+1)) are computed HERE
    (scream computes logn0r from cdistr in a separate helper).
  * both routines leave mu_c/mu_r untouched in the no-condensate branch
    (Fortran intent(out) scalars simply not assigned); callers merge
    with the previous value where that matters (part3, sedimentation).
    Here the masked-out return is 0, matching the zero-initialized
    workspaces at every call site in p3_main.

iparam = 3 (Khairoutdinov-Kogan): the Seifert-Beheng dnu interpolation
branch is dead code and nu stays 0.
"""

import jax.numpy as jnp
from jax.scipy.special import gammaln

from . import constants as c

def _cbrt(x):
    """bfb_cbrt: x**(1/3) via pow (gfortran x**(1._rtype/3._rtype)) —
    NOT libm cbrt, which rounds differently in the last ulp."""
    return x ** (1.0 / 3.0)



def _tgamma(x):
    return jnp.exp(gammaln(x))


def get_cloud_dsd2(qc, nc, rho, context):
    """Cloud DSD parameters (get_cloud_dsd2).

    Returns (nc, mu_c, nu, lamc, cdist, cdist1); nc is updated (floored
    at nsmall and recomputed where the lambda limiters bind).
    """
    qc = jnp.asarray(qc)
    nc = jnp.asarray(nc)
    rho = jnp.asarray(rho)

    gt = (qc >= c.qsmall) & context

    nc_lim = jnp.where(gt, jnp.maximum(nc, c.nsmall), 1.0)
    mu_l = 0.0005714 * (nc_lim * 1.0e-6 * rho) + 0.2714
    mu_l = 1.0 / (mu_l * mu_l) - 1.0
    mu_c = jnp.where(gt, jnp.clip(mu_l, 2.0, 15.0), 0.0)

    qc_safe = jnp.where(gt, qc, 1.0)
    lamc = jnp.where(
        gt,
        _cbrt(c.cons1 * nc_lim * (mu_c + 3.0) * (mu_c + 2.0)
                 * (mu_c + 1.0) / qc_safe),
        0.0)

    lammin = (mu_c + 1.0) * 2.5e4
    lammax = (mu_c + 1.0) * 1.0e6
    lt_min = (lamc < lammin) & gt
    gt_max = (lamc > lammax) & gt
    lamc = jnp.where(lt_min, lammin, jnp.where(gt_max, lammax, lamc))

    either = lt_min | gt_max
    nc_new = jnp.where(
        either,
        6.0 * (lamc * lamc * lamc) * qc_safe
        / (c.pi * c.rho_h2o * (mu_c + 3.0) * (mu_c + 2.0) * (mu_c + 1.0)),
        jnp.where(gt, nc_lim, nc))

    lamc_safe = jnp.where(gt, lamc, 1.0)
    cdist = jnp.where(gt, nc_new * (mu_c + 1.0) / lamc_safe, 0.0)
    cdist1 = jnp.where(gt, nc_new / _tgamma(mu_c + 1.0), 0.0)
    nu = jnp.zeros_like(qc)  # iparam != 1
    return nc_new, mu_c, nu, lamc, cdist, cdist1


def get_rain_dsd2(qr, nr, p3_max_mean_rain_size, context):
    """Rain DSD parameters (get_rain_dsd2, EAM form).

    Returns (nr, mu_r, lamr, cdistr, logn0r); nr is updated (floored at
    nsmall where rain present, recomputed where the limiters bind).
    """
    qr = jnp.asarray(qr)
    nr = jnp.asarray(nr)

    gt = (qr >= c.qsmall) & context

    nr_lim = jnp.maximum(nr, c.nsmall)
    mu_r = jnp.where(gt, c.mu_r_constant, jnp.zeros_like(qr))

    qr_safe = jnp.where(gt, qr, 1.0)
    nr_safe = jnp.where(gt, nr_lim, 1.0)
    lamr = jnp.where(gt,
                     _cbrt(c.cons1 * nr_safe * (mu_r + 3.0) * (mu_r + 2.0)
                              * (mu_r + 1.0) / qr_safe),
                     0.0)

    lammax = (mu_r + 1.0) * 1.0e5
    lammin = (mu_r + 1.0) * (1.0) / p3_max_mean_rain_size
    lt = gt & (lamr < lammin)
    gtm = gt & (lamr > lammax)
    lamr = jnp.where(lt, lammin, jnp.where(gtm, lammax, lamr))

    lamr_safe = jnp.where(gt & (lamr > 0), lamr, 1.0)
    nr_recomp = jnp.exp(3.0 * jnp.log(lamr_safe) + jnp.log(qr_safe)
                        + jnp.log(_tgamma(mu_r + 1.0))
                        - jnp.log(_tgamma(mu_r + 4.0))) / c.cons1
    nr_new = jnp.where(gt, nr_lim, nr)
    nr_new = jnp.where(lt | gtm, nr_recomp, nr_new)

    nr_out_safe = jnp.where(gt & (nr_new > 0), nr_new, 1.0)
    cdistr = jnp.where(gt, nr_new / _tgamma(mu_r + 1.0), 0.0)
    logn0r = jnp.where(gt,
                       jnp.log10(nr_out_safe)
                       + (mu_r + 1.0) * jnp.log10(lamr_safe)
                       - jnp.log10(_tgamma(mu_r + 1.0)), 0.0)
    return nr_new, mu_r, lamr, cdistr, logn0r


def calc_bulk_rho_rime(qi_tot, qi_rim, bi_rim, context):
    """Bulk rime density with limiters (calc_bulkRhoRime).
    Returns (rho_rime, qi_rim, bi_rim)."""
    qi_tot = jnp.asarray(qi_tot)
    qi_rim = jnp.asarray(qi_rim)
    bi_rim = jnp.asarray(bi_rim)

    gt = (bi_rim >= 1.0e-15) & context
    lt = (bi_rim < 1.0e-15) & context

    bi_safe = jnp.where(gt, bi_rim, 1.0)
    rho_rime = jnp.where(gt, qi_rim / bi_safe, 0.0)

    lo = rho_rime < c.rho_rimeMin
    hi = rho_rime > c.rho_rimeMax
    rho_rime = jnp.where(gt & lo, c.rho_rimeMin, rho_rime)
    rho_rime = jnp.where(gt & hi, c.rho_rimeMax, rho_rime)
    adjust = gt & (lo | hi)
    rho_safe = jnp.where(rho_rime > 0, rho_rime, 1.0)
    bi_rim = jnp.where(adjust, qi_rim / rho_safe, bi_rim)

    qi_rim = jnp.where(lt, 0.0, qi_rim)
    bi_rim = jnp.where(lt, 0.0, bi_rim)
    rho_rime = jnp.where(lt, 0.0, rho_rime)

    over = (qi_rim > qi_tot) & (rho_rime > 0) & context
    qi_rim = jnp.where(over, qi_tot, qi_rim)
    bi_rim = jnp.where(over, qi_rim / rho_safe, bi_rim)

    tiny = (qi_rim < c.qsmall) & context
    qi_rim = jnp.where(tiny, 0.0, qi_rim)
    bi_rim = jnp.where(tiny, 0.0, bi_rim)
    return rho_rime, qi_rim, bi_rim
