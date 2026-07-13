"""P3 lookup tables: ice tables (read), rain tables (computed), dnu.

Source: micro_p3.F90 p3_init_a (ice/collection table file read) and
p3_init_b (rain fallspeed/ventilation/evaporation tables). COPIED from
scream_jax/p3/tables.py and adapted:

  * mu_r_constant = 0 (EAM) instead of 1 — the rain tables differ
    numerically from SCREAM's.
  * table version is an argument (4.1.1 is the local file; EAMv3's
    default 4.1.2 is not on local disk — pure input data either way).
  * the ventilation-weight line of p3_init_b uses a SINGLE-precision
    1.e-6 literal (dum5 accumulation), unlike the 1.e-6_rtype in the
    number/mass weights; reproduced via np.float32.
  * expressions kept in the Fortran's 10**(mu*log10(dia)+...) form (not
    simplified to dia**mu) so the port matches pow-for-pow.

The golden archive additionally stores the Fortran-generated tables
(p3_get_tables); replay uses those, and a test asserts this
recomputation matches them (summation-order roundoff only).
"""

import functools

import numpy as np

from . import constants as c

DENSIZE, RIMSIZE, ISIZE = c.densize, c.rimsize, c.isize
ICE_TABLE_SIZE = c.ice_table_size
RCOLLSIZE, COLLECT_TABLE_SIZE = c.rcollsize, c.collect_table_size
VTABLE_DIM0, VTABLE_DIM1 = 300, 10
MU_R_TABLE_DIM = 150


def read_ice_lookup_tables(filename, version="4.1.1"):
    """Parse the ice/collection lookup tables (p3_init_a).

    Returns (ice_table_vals, collect_table_vals) with shapes
    (5, 4, 50, 12) and (5, 4, 50, 30, 2); collection values are log10.
    """
    with open(filename) as f:
        raw = f.read().split()
    assert raw[0] == "VERSION", f"Bad {filename}: expected VERSION header"
    assert raw[1] == version, (
        f"Bad {filename}: expected version {version}, got {raw[1]}")

    tokens = np.array(raw[2:], dtype=np.float64)

    ice = np.empty((DENSIZE, RIMSIZE, ISIZE, ICE_TABLE_SIZE))
    coll = np.empty((DENSIZE, RIMSIZE, ISIZE, RCOLLSIZE, COLLECT_TABLE_SIZE))

    pos = 0
    ice_row = 2 + 15    # 2 int labels + 15 values
    coll_row = 2 + 6    # 2 int labels + 6 values
    # p3_init_a keeps dumk(1..8) then skips one then dumk(9..12):
    # of the 15 values, drop [0], [1] and [10]
    ice_keep = [j for j in range(15) if j > 1 and j != 10]
    for jj in range(DENSIZE):
        for ii in range(RIMSIZE):
            block = tokens[pos:pos + ISIZE * ice_row].reshape(ISIZE, ice_row)
            ice[jj, ii] = block[:, 2:][:, ice_keep]
            pos += ISIZE * ice_row

            block = tokens[pos:pos + ISIZE * RCOLLSIZE * coll_row] \
                .reshape(ISIZE, RCOLLSIZE, coll_row)
            coll[jj, ii] = np.log10(block[:, :, 2:][:, :, [3, 4]])
            pos += ISIZE * RCOLLSIZE * coll_row
    assert pos == tokens.size, "ice lookup table size mismatch"
    return ice, coll


def compute_rain_tables():
    """Recompute the rain fallspeed/ventilation tables (p3_init_b).

    Returns (mu_r_table_vals (150,), vn_table_vals (300,10),
    vm_table_vals (300,10), revap_table_vals (300,10)). All 10 mu_r
    columns are identical since mu_r = mu_r_constant = 0 throughout.
    """
    thrd = c.thrd
    small = 1.0e-30
    mu_r = c.mu_r_constant     # 0.0 in EAM
    dd = 2.0
    # the dum5 weight line uses a single-precision 1.e-6 literal
    w6_single = float(np.float64(np.float32(1.0e-6)))

    mu_r_table = np.full(MU_R_TABLE_DIM, c.mu_r_constant)

    jjs = np.arange(1, VTABLE_DIM0 + 1)
    dm = np.where(jjs <= 20, (jjs * 10.0 - 5.0) * 1e-6,
                  ((jjs - 20) * 30.0 + 195.0) * 1e-6)
    lamr = (mu_r + 1.0) / dm                      # (300,)

    kks = np.arange(1, 10001)
    dia = (kks * dd - dd * 0.5) * 1e-6            # (10000,)
    amg = c.piov6 * 997.0 * dia ** 3 * 1000.0     # mass in [g]
    dia_um = dia * 1e6
    vt = np.where(dia_um <= 134.43, 4.5795e3 * amg ** (2.0 * thrd),
         np.where(dia_um < 1511.64, 4.962e1 * amg ** thrd,
         np.where(dia_um < 3477.84, 1.732e1 * amg ** c.sxth, 9.17)))

    log_dia = np.log10(dia)
    expfac = np.exp(-lamr[:, None] * dia[None, :])         # (300, 10000)
    w_n = 10.0 ** (mu_r * log_dia + 4.0 * mu_r) * dd * 1e-6
    w_m = 10.0 ** ((mu_r + 3.0) * log_dia + 4.0 * mu_r) * dd * 1e-6
    w_v = (vt * dia) ** 0.5 * 10.0 ** ((mu_r + 1.0) * log_dia
                                       + 3.0 * mu_r) * dd * w6_single

    dum1 = (vt * w_n * expfac).sum(axis=1)
    dum2 = np.maximum((w_n * expfac).sum(axis=1), small)
    dum3 = (vt * w_m * expfac).sum(axis=1)
    dum4 = np.maximum((w_m * expfac).sum(axis=1), small)
    dum5 = np.maximum((w_v * expfac).sum(axis=1), small)

    vn_col = dum1 / dum2
    vm_col = dum3 / dum4
    revap_col = 10.0 ** (np.log10(dum5) + (mu_r + 1.0) * np.log10(lamr)
                         - 3.0 * mu_r)

    vn = np.tile(vn_col[:, None], (1, VTABLE_DIM1))
    vm = np.tile(vm_col[:, None], (1, VTABLE_DIM1))
    revap = np.tile(revap_col[:, None], (1, VTABLE_DIM1))
    return mu_r_table, vn, vm, revap


def compute_dnu():
    """Droplet spectral shape parameter array (micro_p3_utils_init)."""
    return c.dnu.copy()


@functools.lru_cache(maxsize=None)
def p3_init(table_dir, version="4.1.1"):
    """Load/compute all P3 lookup tables (p3_init = p3_init_a+p3_init_b).
    Returns a dict of numpy arrays."""
    ice, coll = read_ice_lookup_tables(
        f"{table_dir}/p3_lookup_table_1.dat-v{version}", version)
    mu_r, vn, vm, revap = compute_rain_tables()
    return {
        "ice_table_vals": ice,
        "collect_table_vals": coll,
        "mu_r_table_vals": mu_r,
        "vn_table_vals": vn,
        "vm_table_vals": vm,
        "revap_table_vals": revap,
        "dnu_table_vals": compute_dnu(),
    }
