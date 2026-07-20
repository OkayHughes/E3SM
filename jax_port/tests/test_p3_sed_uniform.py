"""Tests for the "uniform" (fixed-substep, smooth) sedimentation mode
and the p3_soft_masks option (design B1).

Guards, on golden-derived states (p3_218x72_dt1800_5steps.npz):

  1. sed_mode backward compatibility: sed_mode="while" / "scan" are
     bitwise identical to use_while_loop=True / False, and sed_mode
     wins over use_while_loop when both are given;
  2. uniform-mode water conservation (flux-form identity) at the
     default M for each kernel: the total advected mass change equals
     the surface precip mass to ~1e-13 of the column mass;
  3. no-condensate columns are exact no-ops in uniform mode (zero
     fluxes, zero precip contribution) without any hard column mask;
  4. reverse-mode grads through each uniform kernel and through the
     full p3_process_step with sed_mode="uniform" (with and without
     p3_soft_masks, with smooth_width=0.1 all families) are finite and
     NaN-free;
  5. p3_soft_masks=True changes the primal only at the O(qsmall) level
     on the golden state (hard masks off: inactive columns' computed
     contributions are negligible).
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import jax  # noqa: E402
import jax.numpy as jnp  # noqa: E402

import scream_jax  # noqa: E402,F401  (enables x64)
from scream_jax.foundation import constants as c  # noqa: E402
from scream_jax.p3.sedimentation import (  # noqa: E402
    MAX_SEDI_SUBSTEPS_CLOUD,
    MAX_SEDI_SUBSTEPS_ICE,
    MAX_SEDI_SUBSTEPS_RAIN,
)

# reuse the golden fixture and kernel-arg builders of the AD tests
from test_ad_p3 import (  # noqa: E402,F401
    NCOL_GRAD,
    _KERNELS,
    _cloud_args,
    _ice_args,
    _rain_args,
    _sub,
    golden,
)

DEFAULT_M = {
    "cloud": MAX_SEDI_SUBSTEPS_CLOUD,
    "rain": MAX_SEDI_SUBSTEPS_RAIN,
    "ice": MAX_SEDI_SUBSTEPS_ICE,
}
_PRT_ACCUMULATES = {"cloud": False, "rain": True, "ice": True}
GRAD_M = 32  # reduced substep count for the (cheaper) kernel grad tests


def _run_kernel(g, name, **over):
    fn, argfn = _KERNELS[name][:2]
    kw = argfn(g)
    kw.update(over)
    return fn(**kw)


# ---------------------------------------------------------------------------
# 1. sed_mode aliases are bitwise equal to the legacy flag
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("name", ["cloud", "rain", "ice"])
def test_sed_mode_backward_compat(golden, name):
    legacy_while = _run_kernel(golden, name, use_while_loop=True)
    legacy_scan = _run_kernel(golden, name, use_while_loop=False)
    mode_while = _run_kernel(golden, name, sed_mode="while")
    mode_scan = _run_kernel(golden, name, sed_mode="scan")
    # sed_mode wins over a contradicting use_while_loop
    mode_wins = _run_kernel(golden, name, use_while_loop=True,
                            sed_mode="scan")
    for k in legacy_while:
        assert np.array_equal(np.asarray(mode_while[k]),
                              np.asarray(legacy_while[k]),
                              equal_nan=True), f"{name}/while alias: {k}"
        assert np.array_equal(np.asarray(mode_scan[k]),
                              np.asarray(legacy_scan[k]),
                              equal_nan=True), f"{name}/scan alias: {k}"
        assert np.array_equal(np.asarray(mode_wins[k]),
                              np.asarray(legacy_scan[k]),
                              equal_nan=True), f"{name}/sed_mode wins: {k}"


# ---------------------------------------------------------------------------
# 2. uniform-mode water conservation (flux-form identity) at default M
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("name", ["cloud", "rain", "ice"])
def test_uniform_conservation(golden, name):
    fn, argfn, qname, _, _, prt = _KERNELS[name]
    kw = argfn(golden)
    out = fn(**kw, sed_mode="uniform")

    inv_rho = np.asarray(kw["inv_rho"])
    inv_dz = np.asarray(kw["inv_dz"])
    dt = golden["dt"]

    def colmass(q):
        return np.sum(np.asarray(q) / (inv_rho * inv_dz), axis=-1)

    m0 = colmass(kw[qname])
    lost = m0 - colmass(out[qname])
    contrib = np.asarray(out[prt])
    if _PRT_ACCUMULATES[name]:
        contrib = contrib - np.asarray(kw[f"{prt}_in"])
    surf = contrib * c.RHO_H2O * dt

    assert np.asarray(out["converged"]).all()
    assert np.all(np.isfinite(np.asarray(out[qname])))
    # flux-form identity: mass change == surface precip mass, to FP
    # accumulation error (measured <= 2.4e-15 of the column mass)
    np.testing.assert_array_less(
        np.abs(lost - surf), 1e-13 * np.maximum(m0, 1e-30),
        err_msg=f"{name} uniform mode: flux-form mass identity violated")


# ---------------------------------------------------------------------------
# 3. no-condensate columns: exact no-op without a hard column mask
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("name", ["cloud", "rain", "ice"])
def test_uniform_empty_column_exact_noop(golden, name):
    fn, argfn, qname, qincld, cldname, prt = _KERNELS[name]
    kw = argfn(golden, n=4)
    zero = jnp.zeros_like(jnp.asarray(kw[qname]))
    kw[qname] = zero
    kw[qincld] = zero
    out = fn(**kw, sed_mode="uniform", max_substeps=8)
    assert np.all(np.asarray(out[qname]) == 0.0)
    contrib = np.asarray(out[prt])
    if _PRT_ACCUMULATES[name]:
        contrib = contrib - np.asarray(kw[f"{prt}_in"])
    assert np.all(contrib == 0.0), (
        f"{name} uniform mode: empty columns produced surface precip")


# ---------------------------------------------------------------------------
# 4a. reverse grads through each uniform kernel: finite, NaN-free
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("name", ["cloud", "rain", "ice"])
def test_uniform_kernel_reverse_grad_finite(golden, name):
    fn, argfn, qname, qincld, cldname, prt = _KERNELS[name]
    kw = argfn(golden, n=NCOL_GRAD)
    kw["sed_mode"] = "uniform"
    kw["max_substeps"] = GRAD_M
    cld = jnp.asarray(kw[cldname])

    def scalar(q):
        kw2 = dict(kw)
        kw2[qname] = q
        kw2[qincld] = q / cld
        out = fn(**kw2)
        return (jnp.sum(out[qname])
                + jnp.sum(out[prt]) * c.RHO_H2O * golden["dt"])

    grad = np.asarray(jax.grad(scalar)(jnp.asarray(kw[qname])))
    assert int(np.isnan(grad).sum()) == 0, f"{name}: NaN in uniform grad"
    assert int(np.isinf(grad).sum()) == 0, f"{name}: Inf in uniform grad"
    assert np.linalg.norm(grad) > 0.0


# ---------------------------------------------------------------------------
# 4b. full p3_process_step with uniform sed (+/- soft masks, smooth 0.1)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("soft", [False, True])
def test_process_uniform_reverse_grad_finite(golden, soft):
    # _process_scalar hard-codes scan mode, so wrap p3_process_step here
    g, s = golden, golden["s"]
    n = 8
    x0 = jnp.asarray(_sub(s["qc"], n))

    def scalar_u(x):
        from scream_jax.p3.process import p3_process_step
        fields = {k: jnp.asarray(_sub(s[k], n)) for k in
                  ("T_mid", "p_mid", "p_dry_mid", "pseudo_density",
                   "pseudo_density_dry", "cldfrac_tot", "qv", "qc", "nc",
                   "qr", "nr", "qi", "qm", "ni", "bm",
                   "qv_prev_micro_step", "T_prev_micro_step",
                   "nc_nuceat_tend", "ni_activated", "inv_qc_relvar",
                   "precip_liq_surf_mass", "precip_ice_surf_mass")}
        fields["qc"] = x
        out = p3_process_step(
            g["dt"], True, False, True, False, False, False, False, False,
            fields["T_mid"], fields["p_mid"], fields["p_dry_mid"],
            fields["pseudo_density"], fields["pseudo_density_dry"],
            fields["cldfrac_tot"],
            fields["qv"], fields["qc"], fields["nc"], fields["qr"],
            fields["nr"], fields["qi"], fields["qm"], fields["ni"],
            fields["bm"], fields["qv_prev_micro_step"],
            fields["T_prev_micro_step"], fields["nc_nuceat_tend"], None,
            fields["ni_activated"], fields["inv_qc_relvar"],
            fields["precip_liq_surf_mass"], fields["precip_ice_surf_mass"],
            g["tbl"], g["opts"], sed_mode="uniform", p3_soft_masks=soft,
            smooth_width=0.1, smooth_families=None)
        return (jnp.sum(out["qc"] + out["qr"]) + jnp.sum(out["T_mid"])
                + jnp.sum(out["precip_liq_surf_mass"]
                          + out["precip_ice_surf_mass"]))

    grad = np.asarray(jax.grad(scalar_u)(x0))
    n_nan = int(np.isnan(grad).sum())
    assert n_nan == 0, (
        f"uniform sed (soft={soft}): {n_nan}/{grad.size} NaN in grad")
    assert int(np.isinf(grad).sum()) == 0
    assert np.linalg.norm(grad) > 0.0


# ---------------------------------------------------------------------------
# 5. soft masks: primal deviation is O(qsmall) on the golden state
# ---------------------------------------------------------------------------
def test_soft_masks_primal_deviation_small(golden):
    from scream_jax.p3.process import p3_process_step
    g, s = golden, golden["s"]

    def run(soft):
        return p3_process_step(
            g["dt"], True, False, True, False, False, False, False, False,
            s["T_mid"], s["p_mid"], s["p_dry_mid"], s["pseudo_density"],
            s["pseudo_density_dry"], s["cldfrac_tot"],
            s["qv"], s["qc"], s["nc"], s["qr"], s["nr"], s["qi"], s["qm"],
            s["ni"], s["bm"], s["qv_prev_micro_step"],
            s["T_prev_micro_step"], s["nc_nuceat_tend"], None,
            s["ni_activated"], s["inv_qc_relvar"],
            s["precip_liq_surf_mass"], s["precip_ice_surf_mass"],
            g["tbl"], g["opts"], p3_soft_masks=soft)

    hard, soft = run(False), run(True)
    for k in ("T_mid", "qv", "qc", "qr", "qi",
              "precip_liq_surf_mass", "precip_ice_surf_mass"):
        a, b = np.asarray(hard[k]), np.asarray(soft[k])
        scale = max(np.abs(a).max(), 1e-30)
        rel = np.abs(a - b).max() / scale
        assert rel < 1e-9, (
            f"p3_soft_masks primal deviation for '{k}' is {rel:.3e} "
            "field-scale (expected O(qsmall))")
