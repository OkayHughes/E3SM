"""EAMxx adapter for the JAX cld_fraction implementation.

This module is the `py_module_name` target for the CldFraction atmosphere
process (host backend). Its `main` signature — including in-place writes to
the output numpy arrays — must match what
eamxx_cld_fraction_process_interface.cpp::run_impl passes via
py_module_call("main", ...), i.e. the calling convention of the reference
cld_fraction_numpy.py in the same source directory.

Usage in an EAMxx input.yaml (see TEST_HARNESS_DESIGN.md section 5, Tier 2):
    cld_fraction:
      py_module_name: cld_fraction_jax
      py_module_path: <repo>/jax_port/scream_jax/adapters/eamxx
      py_backend: host

No physics lives here — see scream_jax/cld_fraction/main.py.
"""

import sys
from pathlib import Path

import numpy as np

# Make the scream_jax package importable when the embedded interpreter loads
# this file directly off py_module_path (three levels up: adapters/eamxx/..).
_PKG_ROOT = str(Path(__file__).resolve().parents[3])
if _PKG_ROOT not in sys.path:
    sys.path.insert(0, _PKG_ROOT)

from scream_jax.cld_fraction import cld_fraction_main  # noqa: E402


def init():
    """Called from CldFraction::initialize_impl. Nothing to set up (x64 is
    enabled at scream_jax import time); warm-up compilation happens lazily on
    the first main() call."""
    pass


def main(ice_threshold, ice_4out_threshold,
         qi, liq_cld_frac,
         ice_cld_frac, tot_cld_frac,
         ice_cld_frac_4out, tot_cld_frac_4out):
    """In-place shim: numpy views in (zero-copy EAMxx fields), numpy writes out."""
    ice, tot, ice4out, tot4out = cld_fraction_main(
        float(ice_threshold), float(ice_4out_threshold),
        np.asarray(qi), np.asarray(liq_cld_frac))

    ice_cld_frac[...] = np.asarray(ice)
    tot_cld_frac[...] = np.asarray(tot)
    ice_cld_frac_4out[...] = np.asarray(ice4out)
    tot_cld_frac_4out[...] = np.asarray(tot4out)
