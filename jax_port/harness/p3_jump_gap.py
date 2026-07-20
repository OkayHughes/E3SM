#!/usr/bin/env python3
"""P3 smoothing experiment on the real p3_process_step.

Stage A — measure the jump (boundary) term that hard AD gradients miss.
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

Stage B — measure how much of that gap the approximation-by-identity
smoothed gradient closes (p3_process_step(smooth_width=w), see
scream_jax.p3.SMOOTH_FAMILIES for the adopted jump-site families):

  E3(w)  = <grad J_smooth(qc0; w), qc0*v>      (single evaluation, no
           noise — the matched-width single-point estimator of
           harness/smoothing_estimator_demo.py)
  E3'(w) = (1/N_B) sum_i <grad J_smooth(qc_i; w), qc0*v>  (small
           ensemble under the same noise, first N_B CRN members)

for a width sweep, plus a family ablation (cumulative in the task's
a->e order, and per-family singletons) at one representative width.
Closure is reported against the Stage A E2 truth band.

Run: cd jax_port && .venv/bin/python harness/p3_jump_gap.py
(~2 min Stage A + ~10 min Stage B; each (width, families) config
recompiles the jacrev once — widths are static.)
"""
import argparse
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
N = 64            # Stage A ensemble size
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

DIR = qc0 * v                            # perturbation direction in qc


def run(qc, sed_while, smooth_width=0.0, smooth_families=None):
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
        tbl, opts, sed_use_while_loop=sed_while,
        smooth_width=smooth_width, smooth_families=smooth_families)
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


def stage_a(t0):
    """Stage A. Returns (E1 per-sample arrays, E2 per-sample arrays)."""
    print(f"P3 jump-gap experiment: ncol={NCOL}, N={N}, eta={ETA}, dt={DT}")

    # --- zero-noise reference: jvp and plain FD ---
    _, jvp_all = jax.jvp(lambda qc: run(qc, True), (jnp.asarray(qc0),),
                         (jnp.asarray(DIR),))
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
            E1[name][i] = float(jnp.vdot(g, jnp.asarray(DIR)))
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
    return E1, E2


# --------------------------------------------------------------------------
# Stage B
# --------------------------------------------------------------------------
WIDTHS = (0.02, 0.05, 0.1, 0.2, 0.5)
N_B = 16                # smoothed-gradient ensemble size
ABL_WIDTH = 0.1         # representative width for the family ablation
# cumulative adoption in the task's (a)->(e) family order; None == all
CUMULATIVE = (("homog",),
              ("homog", "tmelt"),
              ("homog", "tmelt", "frz"),
              ("homog", "tmelt", "frz", "evap"),
              ("homog", "tmelt", "frz", "evap", "rime"),
              None)
SINGLETONS = (("tmelt",), ("frz",), ("evap",), ("rime",), ("nucl",))


def fam_label(fams):
    return "ALL" if fams is None else "+".join(fams)


def make_grad3(width, fams):
    def f(qc):
        r = run(qc, False, width, fams)
        return jnp.stack([r[name] for name in OBJS])
    return jax.jit(jax.jacrev(f))


def smoothed_dots(width, fams, n_ens, t0, tag):
    """Returns (E3 dict, E3' mean dict, E3' se dict) for one config."""
    g3 = make_grad3(width, fams)
    tic = time.time()
    g0 = np.asarray(g3(jnp.asarray(qc0)))          # compiles here
    e3 = {name: float(np.vdot(g0[j], DIR)) for j, name in enumerate(OBJS)}
    ens = np.zeros((n_ens, len(OBJS)))
    for i in range(n_ens):
        gi = np.asarray(g3(qc_at(0.0, ETA, XIS[i])))
        ens[i] = [np.vdot(gi[j], DIR) for j in range(len(OBJS))]
    e3p = {name: float(ens[:, j].mean()) for j, name in enumerate(OBJS)}
    e3p_se = {name: float(ens[:, j].std() / np.sqrt(n_ens))
              for j, name in enumerate(OBJS)}
    print(f"  [{time.time()-t0:5.0f}s] {tag:28s} "
          f"(compile+evals {time.time()-tic:.0f}s)")
    return e3, e3p, e3p_se


def report_config(name, e1m, e2lo, e2hi, e3, e3p, e3p_se):
    return (f"    {name:7s} E3={e3[name]:+.4e}  "
            f"E3'={e3p[name]:+.4e} +/- {e3p_se[name]:.1e}   "
            f"[E1={e1m[name]:+.4e} | E2 band {e2lo[name]:+.4e}"
            f" .. {e2hi[name]:+.4e}]")


def stage_b(E1, E2, t0):
    e1m = {name: E1[name].mean() for name in OBJS}
    # E2 truth band: min/max over eps of the ensemble means
    e2m = {name: [E2[name][eps].mean() for eps in EPS_LIST] for name in OBJS}
    e2lo = {name: min(e2m[name]) for name in OBJS}
    e2hi = {name: max(e2m[name]) for name in OBJS}

    print(f"\n=== Stage B: smoothed-gradient closure "
          f"(N_B={N_B}, widths {WIDTHS}) ===")
    results = {}

    print("\n-- width sweep, all families --")
    for w in WIDTHS:
        e3, e3p, se = smoothed_dots(w, None, N_B, t0, f"w={w} ALL")
        results[(w, None)] = (e3, e3p, se)

    print("\n-- family ablation at w=%s: cumulative (a->e order) --"
          % ABL_WIDTH)
    for fams in CUMULATIVE:
        key = (ABL_WIDTH, fams)
        if key in results:      # ALL at ABL_WIDTH reused from the sweep
            continue
        e3, e3p, se = smoothed_dots(ABL_WIDTH, fams, N_B, t0,
                                    f"w={ABL_WIDTH} {fam_label(fams)}")
        results[key] = (e3, e3p, se)

    print("\n-- family ablation at w=%s: singletons --" % ABL_WIDTH)
    for fams in SINGLETONS:
        e3, e3p, se = smoothed_dots(ABL_WIDTH, fams, N_B, t0,
                                    f"w={ABL_WIDTH} {fam_label(fams)}")
        results[(ABL_WIDTH, fams)] = (e3, e3p, se)

    # ---------------- closure tables ----------------
    print("\n================ Stage B closure tables ================")
    print("(closure %% = (E3' - E1) / (mid(E2 band) - E1); E1/E2 from "
          "Stage A, N=%d)" % N)
    for name in OBJS:
        mid = 0.5 * (e2lo[name] + e2hi[name])
        gap = mid - e1m[name]
        print(f"\n  objective '{name}':  E1={e1m[name]:+.4e}   "
              f"E2 band [{e2lo[name]:+.4e}, {e2hi[name]:+.4e}]   "
              f"gap={gap:+.3e}")
        print(f"    {'config':22s} {'E3 (no noise)':>14s} "
              f"{'E3p (N=16)':>14s} {'+/-':>8s} {'closure%':>9s} "
              f"{'sign(E2)':>8s}")
        for (w, fams), (e3, e3p, se) in sorted(
                results.items(),
                key=lambda kv: (kv[0][1] is not None, kv[0])):
            clos = 100.0 * (e3p[name] - e1m[name]) / gap if gap != 0 else 0.0
            samesign = "yes" if np.sign(e3p[name]) == np.sign(mid) else "NO"
            print(f"    w={w:<4g} {fam_label(fams):15s} "
                  f"{e3[name]:>+14.4e} {e3p[name]:>+14.4e} "
                  f"{se[name]:>8.1e} {clos:>8.1f}% {samesign:>8s}")
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=("a", "all"), default="all",
                    help="'a' reproduces Stage A only; 'all' adds Stage B")
    args = ap.parse_args()

    t0 = time.time()
    E1, E2 = stage_a(t0)
    if args.stage == "all":
        stage_b(E1, E2, t0)
    print(f"\ntotal {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
