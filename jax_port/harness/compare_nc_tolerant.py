#!/usr/bin/env python3
"""Tolerant NetCDF comparison for JAX-vs-C++ swap tests.

Like eamxx's compare-nc-files, but with the acceptance criterion of
jax_port/TEST_HARNESS_DESIGN.md §7.2: for each variable, at least
(1 - max_frac) of points must match within --tol (relative to the
variable's max magnitude), and NO point may exceed --hard-tol. This
tolerates isolated knife-edge branch flips (which diverge locally under
last-bit input differences) while still catching real regressions.

Usage:
  compare_nc_tolerant.py -s SRC.nc -t TGT.nc --tol 1e-6 --hard-tol 1e-2 \
      --max-frac 1e-3 -v T_mid qv qc ...
"""

import argparse
import sys

import numpy as np
import netCDF4


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-s", "--src", required=True)
    ap.add_argument("-t", "--tgt", required=True)
    ap.add_argument("--tol", type=float, default=1e-6)
    ap.add_argument("--hard-tol", type=float, default=1e-2)
    ap.add_argument("--max-frac", type=float, default=1e-3)
    ap.add_argument("-v", "--vars", nargs="+", required=True)
    args = ap.parse_args()

    src = netCDF4.Dataset(args.src)
    tgt = netCDF4.Dataset(args.tgt)

    failed = False
    for name in args.vars:
        a = np.array(src[name][:], dtype=np.float64)
        b = np.array(tgt[name][:], dtype=np.float64)
        scale = max(np.abs(a).max(), 1e-30)
        rel = np.abs(a - b) / scale
        frac = float((rel > args.tol).mean())
        worst = float(rel.max())
        ok = frac <= args.max_frac and worst <= args.hard_tol
        status = "PASS" if ok else "FAIL"
        print(f"  {status} {name:20s} max_rel={worst:.3e} "
              f"frac>{args.tol:g}: {frac:.2e} (n={a.size})")
        failed |= not ok

    src.close()
    tgt.close()
    print("Comparisons result:", "FAIL" if failed else "PASS")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
