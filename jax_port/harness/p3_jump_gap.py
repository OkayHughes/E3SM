#!/usr/bin/env python3
"""Stage A of the P3 smoothing experiment: measure the jump (boundary)
term that hard AD gradients miss, on the real p3_process_step.

For a scalar objective J and a relative perturbation direction v in qc,
define F(t) = E_xi[ J(qc*(1 + t*v + eta*xi)) ] with elementwise
Gaussian state noise eta. Compare, with common random numbers:

  E1 = (1/N) sum_i <grad J(qc_i), qc*v>   (ensemble of hard AD grads,
       scan mode; unbiased for the KINK part of F'(0) only)
  E2 = (1/N) sum_i [J(qc_i + eps*qc*v) - J(qc_i - eps*qc*v)] / (2 eps)
       (CRN central FD of the noise-averaged primal, while mode;
       unbiased for ALL of F'(0) incl. jump terms, up to O(eps^2))

The gap E2 - E1 estimates the total boundary contribution from P3's
jump branches for this objective/noise. Zero-noise jvp and FD are
reported for reference.

Run: cd jax_port && .venv/bin/python harness/p3_jump_gap.py
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import jax
import jax.numpy as jnp

jax.config.update("jax_enable_x64", True)

REPO = Path(__file__).resolve().parents[2]
GOLDEN = REPO / "jax_port" / "golden" / \
    "physics_suite_218x72_dt1800_2steps.npz"
DATA = REPO.parent / "e3sm-inputdata" / "atm" / "scream"

NCOL = 16
N = 64            # ensemble size
ETA = 0.05        # relative state noise on qc
EPS_LIST = (0.005, 0.01, 0.02)   # FD half-steps (relative units)
DT = 300.0

z = np.load(GOLDEN)
meta = json.loads(str(z["__metadata__"]))
s = {name: np.asarray(z[f"{name}__step0"], dtype=np.float64)[:NCOL]
     for name in meta["fields"]}

from scream_jax.p3 import DEFAULT_OPTS  # noqa: E402
from scream_jax.p3 import tables as p3_tables  # noqa: E402
from scream_jax.p3.process import p3_process_step  # noqa: E402

tbl = p3_tables.p3_init(str(DATA / "tables"))
opts = dict(DEFAULT_OPTS)

qc0 = s["qc"]
scale = qc0  # multiplicative perturbations: qc * (1 + ...)

rng = np.random.default_rng(7)
v = rng.normal(size=qc0.shape)
v /= np.linalg.norm(v * scale)          # unit direction in the qc metric
XIS = rng.normal(size=(N,) + qc0.shape)  # common random numbers


def run(qc, sed_while):
    out = p3_process_step(
        DT, True, True, True, False, False, False, False, False,
        jnp.asarray(s["T_mid"]), jnp.asarray(s["p_mid"]),
        jnp.asarray(s["p_dry_mid"]), jnp.asarray(s["pseudo_density"]),
        jnp.asarray(s["pseudo_density_dry"]),
        jnp.asarray(s["cldfrac_tot"]),
        jnp.asarray(s["qv"]), qc, jnp.asarray(s["nc"]),
        jnp.asarray(s["qr"]), jnp.asarray(s["nr"]),
        jnp.asarray(s["qi"]), jnp.asarray(s["qm"]),
        jnp.asarray(s["ni"]), jnp.asarray(s["bm"]),
        jnp.asarray(s["qv_prev_micro_step"]),
        jnp.asarray(s["T_prev_micro_step"]),
        jnp.asarray(s["nc_nuceat_tend"]), jnp.asarray(s["nccn"]),
        jnp.asarray(s["ni_activated"]), jnp.asarray(s["inv_qc_relvar"]),
        jnp.asarray(s["precip_liq_surf_mass"]),
        jnp.asarray(s["precip_ice_surf_mass"]),
        tbl, opts, sed_use_while_loop=sed_while)
    return {
        "warm": jnp.sum(out["qc"] + out["qr"]),
        "heat": jnp.sum(out["T_mid"]),
        "precip": jnp.sum(out["precip_liq_surf_mass"]
                          + out["precip_ice_surf_mass"]),
    }


OBJS = ("warm", "heat", "precip")


@jax.jit
def primal_all(qc):
    return run(qc, True)          # fast while-mode primal


def make_grad(name):
    return jax.jit(jax.grad(lambda qc: run(qc, False)[name]))


def qc_at(t_v, xi_amp, xi):
    return jnp.asarray(qc0 * (1.0 + t_v * v + xi_amp * xi))


t0 = time.time()
print(f"P3 jump-gap experiment: ncol={NCOL}, N={N}, eta={ETA}, dt={DT}")

# --- zero-noise reference: jvp and plain FD ---
_, jvp_all = jax.jvp(lambda qc: run(qc, True), (jnp.asarray(qc0),),
                     (jnp.asarray(qc0 * v),))
print("\nzero-noise directional derivative (jvp) and plain FD:")
for name in OBJS:
    fds = []
    for eps in EPS_LIST:
        fp = primal_all(qc_at(eps, 0.0, 0.0))[name]
        fm = primal_all(qc_at(-eps, 0.0, 0.0))[name]
        fds.append(float((fp - fm) / (2 * eps)))
    print(f"  {name:7s} jvp={float(jvp_all[name]):+.6e}  "
          f"fd(eps={EPS_LIST})={[f'{f:+.4e}' for f in fds]}")

# --- ensemble estimators with CRN ---
grads = {name: make_grad(name) for name in OBJS}
E1 = {name: np.zeros(N) for name in OBJS}
E2 = {name: {eps: np.zeros(N) for eps in EPS_LIST} for name in OBJS}
for i in range(N):
    xi = XIS[i]
    qc_i = qc_at(0.0, ETA, xi)
    for name in OBJS:
        g = grads[name](qc_i)
        E1[name][i] = float(jnp.vdot(g, jnp.asarray(qc0 * v)))
    for eps in EPS_LIST:
        fp = primal_all(qc_at(eps, ETA, xi))
        fm = primal_all(qc_at(-eps, ETA, xi))
        for name in OBJS:
            E2[name][eps][i] = float((fp[name] - fm[name]) / (2 * eps))
    if i in (0, N // 2):
        print(f"  ... sample {i} done ({time.time()-t0:.0f}s)")

print(f"\nnoise-averaged sensitivity F'(0), N={N} "
      f"(mean +/- MC standard error):")
for name in OBJS:
    e1m, e1s = E1[name].mean(), E1[name].std() / np.sqrt(N)
    print(f"  {name:7s} E1 (hard-AD ensemble) = {e1m:+.6e} "
          f"+/- {e1s:.2e}")
    for eps in EPS_LIST:
        arr = E2[name][eps]
        m, se = arr.mean(), arr.std() / np.sqrt(N)
        gap = m - e1m
        gap_se = np.sqrt(se ** 2 + e1s ** 2)
        sig = abs(gap) / gap_se if gap_se > 0 else np.inf
        print(f"          E2 (CRN FD, eps={eps:5.3f})  = {m:+.6e} "
              f"+/- {se:.2e}   gap={gap:+.3e} ({sig:.1f} sigma)")
print(f"\ntotal {time.time()-t0:.0f}s")
