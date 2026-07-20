#!/usr/bin/env python3
"""Reverse-mode AD profile of the full SCREAM physics suite and
multi-step rollouts (the ML-training path), on golden states.

For each horizon n (suite steps, dt=1800 each; 12-24 steps = the 6-12 h
training window): primal and reverse-grad wall time, peak process RSS,
NaN/Inf counts, and the gradient norm of scalar objectives
(sum T_mid out; sum accumulated precip) w.r.t. T_mid and qc of the
initial state — plus a forward-mode (jvp) tangent norm through the
T_mid-out map as the cheap alternative measurement. Consecutive-horizon
norm ratios give the per-step gradient amplification.

The rollout feeds each step's output state back as the next input;
surface fields are never modified by the suite so they stay fixed, and
SPA/orbital inputs are precomputed host-side per step (time-driven).

Usage: cd jax_port && .venv/bin/python harness/rollout_ad_profile.py \
    [--ncol 8] [--horizons 1,2,4,8,16,24] [--checkpoint process] \
    [--objectives T,precip] [--jit-max 8]

--checkpoint: process | subcycle | none (jax.checkpoint granularity in
  scream_jax.driver._step_ad).
--jit-max: largest horizon for which the grad is jitted end-to-end
  (compile time grows with n; beyond this, run eagerly — inner
  per-process jits are still cached across steps).
"""
import argparse
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import jax
import jax.numpy as jnp

jax.config.update("jax_enable_x64", True)

REPO = Path(__file__).resolve().parents[2]
GOLDEN = REPO / "jax_port" / "golden" / "physics_suite_218x72_dt1800_2steps.npz"
DATA = REPO.parent / "e3sm-inputdata" / "atm" / "scream"


def rss_gb():
    ru = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # darwin reports bytes, linux kilobytes
    return ru / 2**30 if sys.platform == "darwin" else ru / 2**20


def load(ncol):
    z = np.load(GOLDEN)
    meta = json.loads(str(z["__metadata__"]))
    s = {name: np.asarray(z[f"{name}__step0"], dtype=np.float64)[:ncol]
         for name in meta["fields"]}
    import netCDF4
    ds = netCDF4.Dataset(DATA / "init"
                         / "screami_unit_tests_ne2np4L72_20220822.nc")
    geo = {k: np.array(ds[k][:]) for k in ("hyam", "hybm", "lat", "lon",
                                           "area")}
    ds.close()
    return s, geo, float(meta["dt"])


def build(ncol, geo):
    from scream_jax.driver import ScreamPhysics
    from scream_jax.foundation.thermo import calculate_dx_from_area
    return ScreamPhysics(
        DATA, geo["hyam"], geo["hybm"], geo["lat"][:ncol], geo["lon"][:ncol],
        np.asarray(calculate_dx_from_area(geo["area"][:ncol],
                                          geo["lat"][:ncol])),
        mac_mic_subcycles=6, year=2021, spa_col_indices=np.arange(ncol))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ncol", type=int, default=8)
    ap.add_argument("--horizons", default="1,2,4,8,16,24")
    ap.add_argument("--checkpoint", default="process",
                    choices=["process", "subcycle", "none"])
    ap.add_argument("--objectives", default="T,precip")
    # NOTE: XLA compile of even ONE jitted suite step is ~20 min on CPU
    # (measured); eager reverse mode is ~80 s/step with the per-process
    # inner jits cached. Default is therefore fully eager.
    ap.add_argument("--jit-max", type=int, default=0)
    ap.add_argument("--skip-jvp", action="store_true")
    args = ap.parse_args()

    horizons = [int(h) for h in args.horizons.split(",")]
    objectives = args.objectives.split(",")
    n_max = max(horizons)
    s0, geo, dt = load(args.ncol)
    phys = build(args.ncol, geo)
    doy0 = 284 + 45000.0 / 86400.0

    # host-side hoists, one list of per-substep SPA dicts per step
    spa_all = [phys.spa_substep_outputs(dt, doy0 + k * dt / 86400.0,
                                        s0["p_mid"])
               for k in range(n_max)]

    def rollout(T0, qc0, n):
        st = dict(s0)
        st["T_mid"] = T0
        st["qc"] = qc0
        for k in range(n):
            st = phys.step(st, dt, k, doy0 + k * dt / 86400.0,
                           ad_mode=True, spa_outs=spa_all[k],
                           ad_checkpoint=args.checkpoint)
        return st

    def make_obj(name, n):
        if name == "T":
            return lambda T0, qc0: jnp.sum(rollout(T0, qc0, n)["T_mid"])

        def precip(T0, qc0):
            st = rollout(T0, qc0, n)
            return jnp.sum(st["precip_liq_surf_mass"]
                           + st["precip_ice_surf_mass"])
        return precip

    T0 = jnp.asarray(s0["T_mid"])
    qc0 = jnp.asarray(s0["qc"])
    rng = np.random.default_rng(0)
    vT = rng.normal(size=T0.shape)
    vT /= np.linalg.norm(vT)
    vT = jnp.asarray(vT)

    print(f"ncol={args.ncol} nlev={T0.shape[1]} dt={dt:.0f}s "
          f"checkpoint={args.checkpoint} horizons={horizons} "
          f"jit for n<={args.jit_max}", flush=True)
    print(f"{'n':>3} {'obj':>7} {'jit':>4} {'primal_s':>9} "
          f"{'grad1_s':>8} {'grad2_s':>8} {'rss_GB':>7} "
          f"{'|g_T|':>10} {'|g_qc|':>10} {'NaN':>5}", flush=True)
    norms = {o: {} for o in objectives}

    for n in horizons:
        use_jit = n <= args.jit_max
        for obj_name in objectives:
            obj = make_obj(obj_name, n)
            gfun = jax.grad(obj, argnums=(0, 1))
            if use_jit:
                pfun = jax.jit(obj)
                gfun = jax.jit(gfun)
            else:
                pfun = obj

            t0 = time.time()
            p = float(pfun(T0, qc0))
            t_primal = time.time() - t0

            t0 = time.time()
            g = gfun(T0, qc0)
            jax.block_until_ready(g)
            t_grad1 = time.time() - t0
            # second call: warm caches (inner-jit VJPs compiled) — the
            # steady-state training-iteration cost
            t0 = time.time()
            g = gfun(T0, qc0)
            jax.block_until_ready(g)
            t_grad2 = time.time() - t0

            gT, gq = (np.asarray(x) for x in g)
            nan = int(np.isnan(gT).sum() + np.isnan(gq).sum()
                      + np.isinf(gT).sum() + np.isinf(gq).sum())
            nT, nq = np.linalg.norm(gT), np.linalg.norm(gq)
            norms[obj_name][n] = (nT, nq)
            print(f"{n:>3} {obj_name:>7} {str(use_jit):>4} "
                  f"{t_primal:>9.1f} {t_grad1:>8.1f} {t_grad2:>8.2f} "
                  f"{rss_gb():>7.2f} {nT:>10.3e} {nq:>10.3e} "
                  f"{nan:>5}   primal={p:.6e}", flush=True)

        if not args.skip_jvp and "T" in objectives:
            def tout(T):
                return rollout(T, qc0, n)["T_mid"]
            t0 = time.time()
            _, tang = jax.jvp(tout, (T0,), (vT,))
            t_jvp = time.time() - t0
            tang = np.asarray(tang)
            print(f"    jvp(T->T_out) n={n}: |J v| = "
                  f"{np.linalg.norm(tang):.3e} "
                  f"NaN={np.isnan(tang).sum()} ({t_jvp:.1f}s eager)",
                  flush=True)

    print("\ngradient-norm growth vs horizon (per-step amplification "
          "= ratio^(1/dn)):")
    for obj_name in objectives:
        print(f"  objective sum({obj_name}):")
        prev = None
        for n in horizons:
            nT, nq = norms[obj_name][n]
            amp = ""
            if prev is not None:
                dn = n - prev[0]
                if prev[1] > 0:
                    amp = (f"  ratio |g_T| x{nT / prev[1]:.3f} "
                           f"(per-step x{(nT / prev[1]) ** (1 / dn):.3f})")
            print(f"    n={n:>2}: |g_T|={nT:.4e} |g_qc|={nq:.4e}{amp}")
            prev = (n, nT)


if __name__ == "__main__":
    main()
