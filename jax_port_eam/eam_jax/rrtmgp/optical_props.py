"""Optical property operations: delta scaling and increments.

COPIED from ../jax_port/scream_jax/rrtmgp/optical_props.py and adapted
to the EAM FORTRAN kernels
(components/eam/src/physics/rrtmgp/external/rte/kernels/
mo_optical_props_kernels.F90):

  * delta_scale_2str_k (Fortran) has NO `tau > eps` gate — the scaling
    is applied at every point — and guards the two divisions with
    max(eps, 1-wf) / max(eps, 1-f). The C++/EAMxx variant the SCREAM
    port followed gates on tau > eps and divides unguarded; with g = 0
    outside clouds both give tau'=tau, ssa'=ssa, g'=0, but the Fortran
    form is transcribed exactly here.
  * delta_scale on 1-scalar (LW) optics is a no-op in
    mo_optical_props.F90 (delta_scale_1scl) — nothing to port.
  * increment kernels are identical between the two ancestries
    (eps = 3*tiny; max(eps, .) guarded divisions).

Optical properties are dicts: 2-stream {tau, ssa, g}, 1-scalar {tau},
arrays (ncol, nlay, ngpt-or-nbnd).
"""

import jax.numpy as jnp

_EPS = 3.0 * 2.2250738585072014e-308


def delta_scale_2str(op):
    """delta_scale_2str_k (Fortran variant, f = g*g, no tau gate)."""
    tau, ssa, g = op["tau"], op["ssa"], op["g"]
    f = g * g
    wf = ssa * f
    tau_new = (1.0 - wf) * tau
    ssa_new = (ssa - wf) / jnp.maximum(_EPS, 1.0 - wf)
    g_new = (g - f) / jnp.maximum(_EPS, 1.0 - f)
    return {"tau": tau_new, "ssa": ssa_new, "g": g_new}


def _expand_bybnd(arr_bnd, gpt2band):
    """Map a (..., nbnd) array to (..., ngpt) via the band of each gpt."""
    return arr_bnd[..., jnp.asarray(gpt2band)]


def increment_2stream_by_2stream(op1, op2, gpt2band=None):
    """op1 (by g-point) incremented by op2; if op2 is by band, pass
    gpt2band to expand it (inc_2stream_by_2stream_bybnd)."""
    tau1, ssa1, g1 = op1["tau"], op1["ssa"], op1["g"]
    tau2, ssa2, g2 = op2["tau"], op2["ssa"], op2["g"]
    if gpt2band is not None:
        tau2 = _expand_bybnd(tau2, gpt2band)
        ssa2 = _expand_bybnd(ssa2, gpt2band)
        g2 = _expand_bybnd(g2, gpt2band)
    tau12 = tau1 + tau2
    tauscat12 = tau1 * ssa1 + tau2 * ssa2
    g_new = (tau1 * ssa1 * g1 + tau2 * ssa2 * g2) / jnp.maximum(_EPS, tauscat12)
    ssa_new = tauscat12 / jnp.maximum(_EPS, tau12)
    return {"tau": tau12, "ssa": ssa_new, "g": g_new}


def increment_1scalar_by_1scalar(op1, op2, gpt2band=None):
    """op1 tau incremented by op2 tau (by band if gpt2band given)."""
    tau2 = op2["tau"]
    if gpt2band is not None:
        tau2 = _expand_bybnd(tau2, gpt2band)
    return {"tau": op1["tau"] + tau2}
