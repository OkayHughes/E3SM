"""P3 lookup tables: ice tables (read), rain tables (computed), dnu.

Source: components/eamxx/src/physics/p3/impl/p3_init_impl.hpp and the
table declarations in p3_functions.hpp:

    densize=5, rimsize=4, isize=50, ice_table_size=12,
    rcollsize=30, collect_table_size=2
    rain tables: (300, 10) [VTABLE_DIM0 x VTABLE_DIM1], mu_r table: 150

The ice table file (scream/tables/p3_lookup_table_1.dat-v4.1.1) is
whitespace-separated text with a "VERSION 4.1.1" header. The rain tables
are recomputed here exactly as in compute_tables (they can also be read
from the binary .dat8 caches; regeneration matches those to roundoff — a
test asserts it).

All loading is host-side numpy (done once at init); kernels receive plain
arrays.
"""

import functools

import numpy as np

from ..foundation import constants as c

DENSIZE, RIMSIZE, ISIZE = 5, 4, 50
ICE_TABLE_SIZE = 12
RCOLLSIZE, COLLECT_TABLE_SIZE = 30, 2
VTABLE_DIM0, VTABLE_DIM1 = 300, 10
MU_R_TABLE_DIM = 150
P3_VERSION = "4.1.1"


def read_ice_lookup_tables(filename):
    """Parse the ice/collection lookup tables (read_ice_lookup_tables).

    Returns (ice_table_vals, collect_table_vals) with shapes
    (5, 4, 50, 12) and (5, 4, 50, 30, 2); collection values are log10.
    """
    with open(filename) as f:
        raw = f.read().split()
    assert raw[0] == "VERSION", f"Bad {filename}: expected VERSION header"
    assert raw[1] == P3_VERSION, (
        f"Bad {filename}: expected version {P3_VERSION}, got {raw[1]}")

    # Rows are ragged (ice rows: 2 labels + 15 values; collection rows:
    # 2 labels + 6 values), so stream tokens rather than loadtxt.
    tokens = np.array(raw[2:], dtype=np.float64)

    ice = np.empty((DENSIZE, RIMSIZE, ISIZE, ICE_TABLE_SIZE))
    coll = np.empty((DENSIZE, RIMSIZE, ISIZE, RCOLLSIZE, COLLECT_TABLE_SIZE))

    pos = 0
    ice_row = 2 + 15    # 2 int labels + 15 values
    coll_row = 2 + 6    # 2 int labels + 6 values
    # value columns kept: j > 1 and j != 10 -> indices 2..9, 11..14 (12 vals)
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
    """Recompute the rain fallspeed/ventilation tables (compute_tables).

    Returns (mu_r_table_vals (150,), vn_table_vals (300,10),
    vm_table_vals (300,10), revap_table_vals (300,10)).
    Vectorized transcription of the C++ triple loop (mu_r is constant 1
    in table version 4, so all 10 columns are identical by construction).
    """
    thrd = 1.0 / 3.0
    small = 1.0e-30
    mu_r = 1.0
    dd = 2.0

    mu_r_table = np.ones(MU_R_TABLE_DIM)

    # Mean-size axis (jj = 1..300)
    jjs = np.arange(1, VTABLE_DIM0 + 1)
    dm = np.where(jjs <= 20, (jjs * 10 - 5) * 1e-6,
                  ((jjs - 20) * 30 + 195) * 1e-6)
    lamr = (mu_r + 1.0) / dm                      # (300,)

    # PSD bins (kk = 1..10000)
    kks = np.arange(1, 10001)
    dia = (kks * dd - dd * 0.5) * 1e-6            # (10000,)
    amg = c.PIOV6 * 997.0 * dia ** 3 * 1000.0     # mass in [g]
    dia_um = dia * 1e6
    vt = np.where(dia_um <= 134.43, 4.5795e3 * amg ** (2 * thrd),
         np.where(dia_um < 1511.64, 4.962e1 * amg ** thrd,
         np.where(dia_um < 3477.84, 1.732e1 * amg ** c.SXTH, 9.17)))

    log_dia = np.log10(dia)
    expfac = np.exp(-lamr[:, None] * dia[None, :])         # (300, 10000)
    w_n = 10.0 ** (mu_r * log_dia + 4 * mu_r) * dd * 1e-6  # (10000,)
    w_m = 10.0 ** ((mu_r + 3) * log_dia + 4 * mu_r) * dd * 1e-6
    w_v = (vt * dia) ** 0.5 * 10.0 ** ((mu_r + 1) * log_dia + 3 * mu_r) * dd * 1e-6

    dum1 = (vt * w_n * expfac).sum(axis=1)
    dum2 = np.maximum((w_n * expfac).sum(axis=1), small)
    dum3 = (vt * w_m * expfac).sum(axis=1)
    dum4 = np.maximum((w_m * expfac).sum(axis=1), small)
    dum5 = np.maximum((w_v * expfac).sum(axis=1), small)

    vn_col = dum1 / dum2
    vm_col = dum3 / dum4
    revap_col = 10.0 ** (np.log10(dum5) + (mu_r + 1) * np.log10(lamr) - 3 * mu_r)

    vn = np.tile(vn_col[:, None], (1, VTABLE_DIM1))
    vm = np.tile(vm_col[:, None], (1, VTABLE_DIM1))
    revap = np.tile(revap_col[:, None], (1, VTABLE_DIM1))
    return mu_r_table, vn, vm, revap


def read_computed_tables(table_dir):
    """Read the cached binary rain tables (io_impl<true>, double precision)."""
    def rd(name, shape):
        a = np.fromfile(f"{table_dir}/{name}_v2.dat8", dtype=np.float64)
        return a.reshape(shape)
    return (rd("mu_r_table_vals", (MU_R_TABLE_DIM,)),
            rd("vn_table_vals", (VTABLE_DIM0, VTABLE_DIM1)),
            rd("vm_table_vals", (VTABLE_DIM0, VTABLE_DIM1)),
            rd("revap_table_vals", (VTABLE_DIM0, VTABLE_DIM1)))


def compute_dnu():
    """Hard-coded droplet spectral shape parameter array (compute_dnu)."""
    return np.array([0.000, -0.557, -0.430, -0.307, -0.186, -0.067,
                     -0.050, -0.167, -0.282, -0.397, -0.512, -0.626,
                     -0.739, -0.853, -0.966, -0.966])


@functools.lru_cache(maxsize=None)
def p3_init(table_dir):
    """Load/compute all P3 lookup tables (Functions::p3_init).

    table_dir must contain p3_lookup_table_1.dat-v4.1.1; the rain tables
    are recomputed (equivalent to the .dat8 caches to roundoff).
    Returns a dict of numpy arrays.
    """
    ice, coll = read_ice_lookup_tables(
        f"{table_dir}/p3_lookup_table_1.dat-v{P3_VERSION}")
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
