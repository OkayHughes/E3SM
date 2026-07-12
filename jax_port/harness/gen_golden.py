#!/usr/bin/env python3
"""Golden-data generator: run a single EAMxx atmosphere process via pyeamxx
and snapshot every field before/after each timestep into an .npz archive.

The archive is the Tier-0/Tier-1 reference for the JAX ports (see
jax_port/TEST_HARNESS_DESIGN.md §5): state at step n-1 is the input to the
scheme, state at step n is the expected output.

Run INSIDE the scream-dev container, after building with
EAMXX_ENABLE_PYSCREAM=ON (see jax_port/dev/README.md):

    docker exec -w /work/E3SM/jax_port/harness scream-dev \
        python3 gen_golden.py p3 --steps 5 -o /work/E3SM/jax_port/golden/p3.npz

Notes:
- Required fields that are absent from the IC file are filled with the
  per-process constants below (mirroring the single-process input.yaml
  initial_conditions overrides) or zero, and recorded in the archive's
  metadata. Golden data is self-consistent regardless: the archive stores
  the inputs actually used.
- Arrays are stored as '<field>__step<n>' for n = 0 (post-init) .. nsteps.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

BUILD_DIR = "/work/E3SM/components/eamxx/ctest-build/copilot-testing/full_debug"
IC_FILE = "/work/e3sm-inputdata/atm/scream/init/screami_unit_tests_ne2np4L72_20220822.nc"
NCOLS, NLEVS = 218, 72
T0 = "2021-10-12-45000"
DT = 1800

# Process configs: factory params (from the single-process input.yaml files)
# and constant fills for fields missing from the IC file.
PROC_CONFIGS = {
    "p3": {
        "params": {
            "max_total_ni": 740.0e3,
            "do_prescribed_ccn": False,
        },
        "ic_fill": {
            "precip_liq_surf_mass": 0.0,
            "precip_ice_surf_mass": 0.0,
        },
    },
    "shoc": {
        "params": {
            # SHOC requires an internal dt <= 300 s; with the default 1800 s
            # driver step this means 6 subcycles (as in the AD tests).
            "number_of_subcycles": 6,
            "lambda_low": 0.001,
            "lambda_high": 0.08,
            "lambda_slope": 2.65,
            "lambda_thresh": 0.02,
            "thl2tune": 1.0,
            "qw2tune": 1.0,
            "qwthl2tune": 1.0,
            "w2tune": 1.0,
            "length_fac": 0.5,
            "c_diag_3rd_mom": 7.0,
            "coeff_kh": 0.1,
            "coeff_km": 0.1,
            "shoc_1p5tke": False,
        },
        "ic_fill": {
            "surf_sens_flux": 0.0,
            "surf_evap": 0.0,
            "phis": 0.0,  # no topography file in this harness; recorded in metadata
        },
    },
}


def snapshot(proc, names):
    out = {}
    for name in names:
        f = proc.get_field(name)
        f.sync_to_host()
        out[name] = np.array(f.get(), copy=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("process", choices=sorted(PROC_CONFIGS))
    ap.add_argument("--steps", type=int, default=5)
    ap.add_argument("--dt", type=int, default=DT)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--build-dir", default=BUILD_DIR)
    ap.add_argument("--ic-file", default=IC_FILE)
    args = ap.parse_args()

    # The compiled ext lives in <build>/src/python/libpyeamxx (imported as
    # libpyeamxx.pyeamxx_ext); the pure-python `pyeamxx` wrapper lives in the
    # source tree.
    sys.path.append(str(Path(args.build_dir) / "src" / "python"))
    sys.path.append("/work/E3SM/components/eamxx/src/python")

    import mpi4py
    mpi4py.rc.initialize = False
    mpi4py.rc.finalize = False
    from mpi4py import MPI
    import pyeamxx

    MPI.Init()
    pyeamxx.init()
    try:
        run(args, pyeamxx)
    finally:
        pyeamxx.finalize()
        MPI.Finalize()


def run(args, pyeamxx):
    cfg = PROC_CONFIGS[args.process]

    pyeamxx.create_grids_manager(NCOLS, NLEVS, args.ic_file)
    proc = pyeamxx.AtmProc(dict(cfg["params"]), args.process)

    missing = list(proc.read_ic(args.ic_file))
    filled, zero_filled = {}, []
    for name in missing:
        f = proc.get_field(name)
        arr = f.get()
        val = cfg["ic_fill"].get(name, 0.0)
        arr[...] = val
        f.sync_to_dev()
        if name in cfg["ic_fill"]:
            filled[name] = val
        else:
            zero_filled.append(name)
    if zero_filled:
        print(f"WARNING: zero-filled fields missing from IC file: {zero_filled}")

    proc.initialize(T0)

    names = list(proc.list_all_fields())
    arrays = {}
    for name, arr in snapshot(proc, names).items():
        arrays[f"{name}__step0"] = arr

    for n in range(1, args.steps + 1):
        proc.run(args.dt)
        for name, arr in snapshot(proc, names).items():
            arrays[f"{name}__step{n}"] = arr
        print(f"step {n}/{args.steps} captured")

    meta = {
        "process": args.process,
        "params": cfg["params"],
        "ncols": NCOLS,
        "nlevs": NLEVS,
        "dt": args.dt,
        "t0": T0,
        "nsteps": args.steps,
        "ic_file": args.ic_file,
        "ic_fill": filled,
        "zero_filled": zero_filled,
        "fields": names,
        "required_fields": list(proc.list_required_fields()),
        "computed_fields": list(proc.list_computed_fields()),
    }
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out, __metadata__=json.dumps(meta), **arrays)
    print(f"wrote {out} ({len(arrays)} arrays, {len(names)} fields x {args.steps + 1} steps)")


if __name__ == "__main__":
    main()
