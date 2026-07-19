#!/usr/bin/env python3
"""Numerical check of the smoothed-gradient-as-estimator question:

For a jump f(x) = 1{x > tau} and a kink k(x) = max(x - tau, 0) under
Gaussian input noise eta, compare three estimators of
d/dx E_xi[f(x + xi)] (the noise-averaged sensitivity, whose truth is
known analytically):

  (a) ensemble mean of HARD AD gradients  (pathwise, w = 0)
  (b) ensemble mean of SMOOTHED AD gradients at width w  (N samples)
  (c) single-point smoothed gradient at logistic-matched width
      (zero-variance analytic 'Rao-Blackwell' along the switch axis)

Truth: jump -> phi_eta(tau - x); kink -> Phi((x - tau)/eta).
Run: cd jax_port && .venv/bin/python harness/smoothing_estimator_demo.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import jax
import jax.numpy as jnp

jax.config.update("jax_enable_x64", True)
from scream_jax.foundation import smoothing  # noqa: E402

TAU = 0.0
ETA = 0.1          # input noise scale
X0 = TAU + 0.5 * ETA   # evaluate near the switch (worst case)
N = 32             # small ensemble
REPS = 2000        # repetitions to measure estimator mean/std

rng = np.random.default_rng(0)


def grad_step(w):
    return jax.vmap(jax.grad(lambda x: smoothing.step(x - TAU, w)))


def grad_kink_hard():
    return jax.vmap(jax.grad(lambda x: jnp.maximum(x - TAU, 0.0)))


truth_jump = np.exp(-0.5 * ((TAU - X0) / ETA) ** 2) / (
    ETA * np.sqrt(2 * np.pi))
truth_kink = 0.5 * (1 + np.vectorize(__import__("math").erf)(
    (X0 - TAU) / (ETA * np.sqrt(2))))

print(f"x0 = tau + 0.5*eta, eta = {ETA}, N = {N}, reps = {REPS}")
print(f"truth d/dx E[jump] = {truth_jump:.4f}   "
      f"truth d/dx E[kink] = {float(truth_kink):.4f}\n")

# --- kink: hard pathwise estimator (the 'trivial subgradient case') ---
g_kink = grad_kink_hard()
est = np.array([np.mean(np.asarray(
    g_kink(jnp.asarray(X0 + ETA * rng.normal(size=N)))))
    for _ in range(REPS)])
print(f"KINK  hard-AD ensemble: mean {est.mean():.4f} "
      f"(bias {est.mean()-truth_kink:+.4f})  std {est.std():.4f}")

# --- jump: estimators across widths ---
print("\nJUMP  (hard-AD ensemble mean is identically 0 -> bias "
      f"{-truth_jump:.4f}, std 0)")
print(f"{'w':>10} {'mean':>9} {'bias':>9} {'std':>8} {'rmse':>8}")
for w in (0.4, 0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625):
    g = grad_step(w)
    est = np.array([np.mean(np.asarray(
        g(jnp.asarray(X0 + ETA * rng.normal(size=N)))))
        for _ in range(REPS)])
    bias = est.mean() - truth_jump
    rmse = np.sqrt(bias ** 2 + est.var())
    print(f"{w:>10.5f} {est.mean():>9.4f} {bias:>+9.4f} "
          f"{est.std():>8.4f} {rmse:>8.4f}")

# --- (c) zero-sample analytic: logistic kernel matched to the noise ---
w_match = np.sqrt(3) * ETA / np.pi   # logistic std = pi*w/sqrt(3) = eta
g_match = float(grad_step(w_match)(jnp.asarray([X0]))[0])
print(f"\nmatched-width single evaluation (w = sqrt(3)*eta/pi = "
      f"{w_match:.4f}):\n  value {g_match:.4f}  "
      f"(vs truth {truth_jump:.4f}, rel err "
      f"{abs(g_match-truth_jump)/truth_jump:.2%}, variance 0)")
