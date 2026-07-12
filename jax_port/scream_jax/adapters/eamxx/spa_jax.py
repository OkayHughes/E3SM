"""EAMxx adapter for the JAX SPA implementation.

`py_module_name` target for the spa atmosphere process (host backend).
The `main` signature must match the py_module_call in the
EAMXX_HAS_PYTHON branch of eamxx_spa_process_interface.cpp on this
branch. No physics lives here — see scream_jax/spa/process.py.
"""

import sys
from pathlib import Path

import numpy as np

_PKG_ROOT = str(Path(__file__).resolve().parents[3])
if _PKG_ROOT not in sys.path:
    sys.path.insert(0, _PKG_ROOT)

from scream_jax.spa.process import load_spa_data, spa_process_step  # noqa: E402

_DATA = None


def init(spa_data_file):
    global _DATA
    _DATA = load_spa_data(str(spa_data_file))


def main(doy_end_of_step, p_mid, nccn, aero_g_sw, aero_ssa_sw,
         aero_tau_sw, aero_tau_lw):
    assert _DATA is not None, "spa_jax.init() was not called"
    out = spa_process_step(_DATA, float(doy_end_of_step), np.asarray(p_mid))
    nccn[...] = out["nccn"]
    aero_g_sw[...] = out["aero_g_sw"]
    aero_ssa_sw[...] = out["aero_ssa_sw"]
    aero_tau_sw[...] = out["aero_tau_sw"]
    aero_tau_lw[...] = out["aero_tau_lw"]
