"""P3 lookup-table index computation and interpolation.

Source: micro_p3.F90 access_lookup_table[_coll],
find_lookupTable_indices_1a/1b/3. COPIED from
scream_jax/p3/table_lookups.py; adaptations against the EAM Fortran:

  * lookup_rain (find_lookupTable_indices_1b): activity test is
    qr >= qsmall (scream_jax used qr > qsmall).

Index conventions as in the copy: ice-table indices 0-based, rain
Table3 indices kept 1-based with apply_table3 subtracting 1.
"""

import jax.numpy as jnp

from . import constants as c

ISIZE, RIMSIZE, DENSIZE, RCOLLSIZE = c.isize, c.rimsize, c.densize, c.rcollsize

def _cbrt(x):
    """bfb_cbrt: x**(1/3) via pow (gfortran x**(1._rtype/3._rtype)) —
    NOT libm cbrt, which rounds differently in the last ulp."""
    return x ** (1.0 / 3.0)



def lookup_table3(mu_r, lamr, context):
    """Rain-table index location (find_lookupTable_indices_3)."""
    mu_r = jnp.asarray(mu_r)
    lamr_safe = jnp.where(context & (jnp.asarray(lamr) != 0),
                          jnp.asarray(lamr), 1.0)
    dum1 = (mu_r + 1.0) / lamr_safe

    small = dum1 <= 195.0e-6
    # small-size branch (inv_dum3 = 0.1)
    r_small = jnp.clip((dum1 * 1.0e6 + 5.0) * 0.1, 1.0, 20.0)
    i_small = jnp.clip(r_small.astype(jnp.int32), 1, 20)
    # large-size branch (inv_dum3 = thrd*0.1)
    r_large = jnp.clip((dum1 * 1.0e6 - 195.0) * (c.thrd * 0.1) + 20.0,
                       20.0, 300.0)
    i_large = jnp.clip(r_large.astype(jnp.int32), 20, 299)

    rdumii = jnp.where(small, r_small, r_large)
    dumii = jnp.where(small, i_small, i_large)
    dumii = jnp.where(context, dumii, 1)

    rdumjj = jnp.clip(mu_r + 1.0, 1.0, 10.0)
    dumjj = jnp.clip(rdumjj.astype(jnp.int32), 1, 9)
    rdumjj = jnp.where(context, rdumjj, 0.0)
    dumjj = jnp.where(context, dumjj, 1)
    return {"rdumii": rdumii, "dumii": dumii, "rdumjj": rdumjj, "dumjj": dumjj}


def apply_table3(table, tab3):
    """Bilinear interpolation in a (300, 10) rain table. Indices in tab3
    are 1-based (Fortran)."""
    table = jnp.asarray(table)
    ii = tab3["dumii"]
    jj = tab3["dumjj"]
    w_i = tab3["rdumii"] - ii
    w_j = tab3["rdumjj"] - jj

    t00 = table[ii - 1, jj - 1]
    t10 = table[ii, jj - 1]
    t01 = table[ii - 1, jj]
    t11 = table[ii, jj]
    dum1 = t00 + w_i * (t10 - t00)
    dum2 = t01 + w_i * (t11 - t01)
    return dum1 + w_j * (dum2 - dum1)


def lookup_ice(qi, ni, qm, rhop, context):
    """Ice-table index location (find_lookupTable_indices_1a).
    0-based indices."""
    qi = jnp.asarray(qi)
    ni_safe = jnp.where(context & (jnp.asarray(ni) > 0), jnp.asarray(ni), 1.0)
    qi_safe = jnp.where(context & (qi > 0), qi, 1.0e-18)
    qm_safe = jnp.where(context, jnp.asarray(qm), 0.0)
    rhop_safe = jnp.where(context, jnp.asarray(rhop), 400.0)

    dum1 = (jnp.log10(qi_safe / ni_safe) + 18.0) \
        * c.lookup_table_1a_dum1_c - 10.0
    dumi = dum1.astype(jnp.int32)
    dum1 = jnp.clip(dum1, 1.0, float(ISIZE))
    dumi = jnp.clip(dumi, 1, ISIZE - 1)

    dum4 = (qm_safe / qi_safe) * 3.0 + 1.0
    dumii = dum4.astype(jnp.int32)
    dum4 = jnp.clip(dum4, 1.0, float(RIMSIZE))
    dumii = jnp.clip(dumii, 1, RIMSIZE - 1)

    dum5 = jnp.where(rhop_safe <= 650.0,
                     (rhop_safe - 50.0) * 0.005 + 1.0,
                     (rhop_safe - 650.0) * 0.004 + 4.0)
    dumjj = dum5.astype(jnp.int32)
    dum5 = jnp.clip(dum5, 1.0, float(DENSIZE))
    dumjj = jnp.clip(dumjj, 1, DENSIZE - 1)

    return {"dum1": dum1, "dumi": dumi - 1,
            "dum4": dum4, "dumii": dumii - 1,
            "dum5": dum5, "dumjj": dumjj - 1}


def lookup_rain(qr, nr, context):
    """Rain index for the ice-rain collection table
    (find_lookupTable_indices_1b). 0-based dumj.
    EAM activity test: qr >= qsmall (and nr > 0)."""
    qr = jnp.asarray(qr)
    nr = jnp.asarray(nr)
    gt_small = (qr >= c.qsmall) & (nr > 0.0) & context

    nr_safe = jnp.where(gt_small, nr, 1.0)
    qr_safe = jnp.where(gt_small, qr, 1.0)
    dumlr = jnp.where(gt_small,
                      _cbrt(qr_safe / (c.pi * c.rho_h2o * nr_safe)), 1.0)

    dum3 = (jnp.log10(1.0 * dumlr) + 5.0) * 10.70415
    dumj = dum3.astype(jnp.int32)
    dum3 = jnp.clip(dum3, 1.0, float(RCOLLSIZE))
    dumj = jnp.clip(dumj, 1, RCOLLSIZE - 1)

    dumj = jnp.where(gt_small, dumj, 1)
    dum3 = jnp.where(gt_small, dum3, 1.0)
    return {"dum3": dum3, "dumj": dumj - 1}


def apply_table_ice(idx: int, ice_table_vals, ti):
    """Trilinear interpolation in the (5,4,50,12) ice table
    (access_lookup_table) for quantity column idx (0-based; Fortran
    index = idx+1)."""
    t = jnp.asarray(ice_table_vals)[..., idx]
    i, ii, jj = ti["dumi"], ti["dumii"], ti["dumjj"]
    w1 = ti["dum1"] - i - 1   # (dum1 - dumi_1based)
    w4 = ti["dum4"] - ii - 1
    w5 = ti["dum5"] - jj - 1

    def at(djj, dii, di):
        return t[jj + djj, ii + dii, i + di]

    iproc1 = at(0, 0, 0) + w1 * (at(0, 0, 1) - at(0, 0, 0))
    gproc1 = at(0, 1, 0) + w1 * (at(0, 1, 1) - at(0, 1, 0))
    tmp1 = iproc1 + w4 * (gproc1 - iproc1)

    iproc1 = at(1, 0, 0) + w1 * (at(1, 0, 1) - at(1, 0, 0))
    gproc1 = at(1, 1, 0) + w1 * (at(1, 1, 1) - at(1, 1, 0))
    tmp2 = iproc1 + w4 * (gproc1 - iproc1)

    return tmp1 + w5 * (tmp2 - tmp1)


def apply_table_coll(idx: int, collect_table_vals, ti, tr):
    """Quadrilinear interpolation in the (5,4,50,30,2) collection table
    (access_lookup_table_coll)."""
    t = jnp.asarray(collect_table_vals)[..., idx]
    i, ii, jj = ti["dumi"], ti["dumii"], ti["dumjj"]
    j = tr["dumj"]
    w1 = ti["dum1"] - i - 1
    w3 = tr["dum3"] - j - 1
    w4 = ti["dum4"] - ii - 1
    w5 = ti["dum5"] - jj - 1

    def at(djj, dii, di, dj):
        return t[jj + djj, ii + dii, i + di, j + dj]

    def plane(djj, dii):
        d1 = at(djj, dii, 0, 0) + w1 * (at(djj, dii, 1, 0) - at(djj, dii, 0, 0))
        d2 = at(djj, dii, 0, 1) + w1 * (at(djj, dii, 1, 1) - at(djj, dii, 0, 1))
        return d1 + w3 * (d2 - d1)

    tmp1 = plane(0, 0) + w4 * (plane(0, 1) - plane(0, 0))
    tmp2 = plane(1, 0) + w4 * (plane(1, 1) - plane(1, 0))
    return tmp1 + w5 * (tmp2 - tmp1)
