"""AD hardening tests for P3 sedimentation and the full P3 process step.

The CFL substep loops in scream_jax/p3/sedimentation.py run, by default,
as fixed-length masked lax.scans (reverse-mode differentiable) instead of
lax.while_loop. These tests check, on golden-derived states
(p3_218x72_dt1800_5steps.npz, 218 columns, dt = 1800 s):

  1. primal equivalence: the scan output is bitwise identical to the
     original lax.while_loop implementation (still available behind
     use_while_loop=True) for every output of every kernel;
  2. trip-count safety: every column spends its full dt within the static
     substep bounds (the per-column "converged" flag, i.e. final
     dt_left <= tol), both kernel-level and through p3_main;
  3. reverse-mode grads of each kernel run, are finite and NaN-free;
  4. the dot-product identity <grad f, v> == jvp(f; v) holds per kernel;
  5. reverse-mode grad of the full p3_process_step runs NaN-free
     (d(sum qc+qr)/d(qc) and d(sum qi)/d(qi)).
"""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import jax  # noqa: E402
import jax.numpy as jnp  # noqa: E402

import scream_jax  # noqa: E402,F401  (enables x64)
from scream_jax.foundation import constants as c  # noqa: E402
from scream_jax.foundation.thermo import (  # noqa: E402
    calculate_drymmr_from_wetmmr_dp_based,
    calculate_dz,
    calculate_theta_from_T,
    exner_function,
)
from scream_jax.p3 import DEFAULT_OPTS, tables  # noqa: E402
from scream_jax.p3.main import p3_main  # noqa: E402
from scream_jax.p3.main_part1 import p3_main_part1  # noqa: E402
from scream_jax.p3.process import MINCLD, p3_process_step  # noqa: E402
from scream_jax.p3.sedimentation import (  # noqa: E402
    cloud_sedimentation,
    ice_sedimentation,
    rain_sedimentation,
)

REPO = Path(__file__).resolve().parents[2]
GOLDEN = REPO / "jax_port" / "golden" / "p3_218x72_dt1800_5steps.npz"
TDIR = REPO.parent / "e3sm-inputdata" / "atm" / "scream" / "tables"

NCOL_GRAD = 16  # column subset for the (memory-heavier) reverse-mode tests


@pytest.fixture(scope="module")
def golden():
    """Golden step-0 state (218 cols), tables, opts, and the p3_main /
    sedimentation-kernel inputs derived from it exactly as
    process.py's preamble + p3_main_part1 derive them."""
    if not GOLDEN.exists():
        pytest.skip(f"golden archive not found: {GOLDEN}")
    if not (TDIR / f"p3_lookup_table_1.dat-v{tables.P3_VERSION}").exists():
        pytest.skip("P3 tables not available")
    z = np.load(GOLDEN)
    meta = json.loads(str(z["__metadata__"]))
    dt = float(meta["dt"])
    tbl = tables.p3_init(str(TDIR))
    opts = dict(DEFAULT_OPTS)
    opts["max_total_ni"] = meta["params"]["max_total_ni"]

    s = {k.rsplit("__", 1)[0]: np.asarray(z[k], dtype=np.float64)
         for k in z.files if k.endswith("__step0")}

    # ---- process.py preamble (wet -> dry state on the p3_main grid) ----
    dz = calculate_dz(s["pseudo_density"], s["p_mid"], s["T_mid"], s["qv"])
    dry = {k: calculate_drymmr_from_wetmmr_dp_based(
        s[k], s["pseudo_density"], s["pseudo_density_dry"])
        for k in ("qc", "nc", "qr", "nr", "qi", "ni", "qm", "bm", "qv")}
    inv_exner = 1.0 / exner_function(s["p_mid"])
    th_atm = calculate_theta_from_T(s["T_mid"], s["p_mid"])
    cld_frac_l = np.maximum(s["cldfrac_tot"], MINCLD)
    cld_frac_i = cld_frac_l
    cld_frac_r = cld_frac_l.copy()
    cld_frac_r[..., 1:] = np.maximum(cld_frac_r[..., 1:],
                                     s["cldfrac_tot"][..., :-1])

    # ---- p3_main init + part1 (rho, fall-speed factors, in-cloud q) ----
    exner = 1.0 / inv_exner
    qv0 = np.maximum(np.asarray(dry["qv"]), 0.0)
    T_atm = np.asarray(th_atm) * np.asarray(exner)
    zcol = np.zeros_like(T_atm)
    st = p3_main_part1(
        True, False, dt,
        s["p_dry_mid"], s["pseudo_density_dry"], dz,
        s["nc_nuceat_tend"], zcol, inv_exner, exner,
        1.0 / cld_frac_l, 1.0 / cld_frac_i, 1.0 / cld_frac_r,
        T_atm, qv0, th_atm,
        dry["qc"], dry["nc"], dry["qr"], dry["nr"],
        dry["qi"], dry["ni"], dry["qm"], dry["bm"], opts)
    st = {k: np.asarray(v) for k, v in st.items()}

    return dict(s=s, meta=meta, dt=dt, tbl=tbl, opts=opts, dz=dz,
                inv_dz=1.0 / np.asarray(dz), dry=dry, th_atm=th_atm,
                inv_exner=np.asarray(inv_exner),
                cld_frac_l=cld_frac_l, cld_frac_i=cld_frac_i,
                cld_frac_r=cld_frac_r, st=st, zcol=zcol)


def _sub(a, n):
    return np.asarray(a)[:n] if np.asarray(a).ndim else np.asarray(a)


def _cloud_args(g, n=None, use_while_loop=False):
    st, dt = g["st"], g["dt"]
    f = (lambda a: _sub(a, n)) if n else np.asarray
    zc = f(g["zcol"])
    return dict(
        qc_incld=f(st["qc_incld"]), rho=f(st["rho"]),
        inv_rho=f(st["inv_rho"]), cld_frac_l=f(g["cld_frac_l"]),
        acn=f(st["acn"]), inv_dz=f(g["inv_dz"]), dt=dt, inv_dt=1.0 / dt,
        do_predict_nc=True, qc=f(st["qc"]), nc=f(st["nc"]),
        nc_incld=f(st["nc_incld"]), mu_c=zc, lamc=zc,
        qc_tend_in=f(st["qc"]), nc_tend_in=f(st["nc"]),
        precip_liq_surf_in=np.zeros(zc.shape[:-1]),
        use_while_loop=use_while_loop)


def _rain_args(g, n=None, use_while_loop=False):
    st, dt, tbl = g["st"], g["dt"], g["tbl"]
    f = (lambda a: _sub(a, n)) if n else np.asarray
    zc = f(g["zcol"])
    return dict(
        rho=f(st["rho"]), inv_rho=f(st["inv_rho"]),
        rhofacr=f(st["rhofacr"]), cld_frac_r=f(g["cld_frac_r"]),
        inv_dz=f(g["inv_dz"]), qr_incld=f(st["qr_incld"]),
        vn_table_vals=tbl["vn_table_vals"],
        vm_table_vals=tbl["vm_table_vals"], dt=dt, inv_dt=1.0 / dt,
        qr=f(st["qr"]), nr=f(st["nr"]), nr_incld=f(st["nr_incld"]),
        mu_r=zc, lamr=zc,
        precip_liq_flux_in=np.zeros(zc.shape[:-1] + (zc.shape[-1] + 1,)),
        qr_tend_in=f(st["qr"]), nr_tend_in=f(st["nr"]),
        precip_liq_surf_in=np.zeros(zc.shape[:-1]), opts=g["opts"],
        use_while_loop=use_while_loop)


def _ice_args(g, n=None, use_while_loop=False):
    st, dt, tbl = g["st"], g["dt"], g["tbl"]
    f = (lambda a: _sub(a, n)) if n else np.asarray
    zc = f(g["zcol"])
    return dict(
        rho=f(st["rho"]), inv_rho=f(st["inv_rho"]),
        rhofaci=f(st["rhofaci"]), cld_frac_i=f(g["cld_frac_i"]),
        inv_dz=f(g["inv_dz"]), dt=dt, inv_dt=1.0 / dt,
        qi=f(st["qi"]), qi_incld=f(st["qi_incld"]),
        ni=f(st["ni"]), ni_incld=f(st["ni_incld"]),
        qm=f(st["qm"]), qm_incld=f(st["qm_incld"]),
        bm=f(st["bm"]), bm_incld=f(st["bm_incld"]),
        ice_table_vals=tbl["ice_table_vals"],
        qi_tend_in=f(st["qi"]), ni_tend_in=f(st["ni"]),
        precip_ice_surf_in=np.zeros(zc.shape[:-1]), opts=g["opts"],
        use_while_loop=use_while_loop)


_KERNELS = {
    "cloud": (cloud_sedimentation, _cloud_args, "qc", "qc_incld",
              "cld_frac_l", "precip_liq_surf"),
    "rain": (rain_sedimentation, _rain_args, "qr", "qr_incld",
             "cld_frac_r", "precip_liq_surf"),
    "ice": (ice_sedimentation, _ice_args, "qi", "qi_incld",
            "cld_frac_i", "precip_ice_surf"),
}


def _scalar_fn(g, name):
    """Scalar loss for kernel `name`: column mass out + surface precip,
    as a function of the grid-mean mixing ratio (in-cloud value threaded
    consistently as q/cld_frac)."""
    fn, argfn, qname, qincld, cldname, prt = _KERNELS[name]
    kw = argfn(g, n=NCOL_GRAD)
    cld = jnp.asarray(kw[cldname])

    def scalar(q):
        kw2 = dict(kw)
        kw2[qname] = q
        kw2[qincld] = q / cld
        out = fn(**kw2)
        return (jnp.sum(out[qname])
                + jnp.sum(out[prt]) * c.RHO_H2O * g["dt"])

    x0 = jnp.asarray(kw[qname])
    return scalar, x0


# ---------------------------------------------------------------------------
# (v) primal equivalence: masked scan == original while_loop, bitwise
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("name", ["cloud", "rain", "ice"])
def test_primal_equivalence_scan_vs_while(golden, name):
    fn, argfn = _KERNELS[name][:2]
    a = fn(**argfn(golden))
    b = fn(**argfn(golden, use_while_loop=True))
    assert set(a) == set(b)
    for k in a:
        if k == "converged":  # while variant: trivially all True
            continue
        av, bv = np.asarray(a[k]), np.asarray(b[k])
        assert np.array_equal(av, bv, equal_nan=True), (
            f"{name} sedimentation: scan != while for output '{k}' "
            f"(max abs diff {np.nanmax(np.abs(av - bv)):.3e})")


# ---------------------------------------------------------------------------
# (iv) trip-count safety: dt fully spent within the static substep bounds
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("name", ["cloud", "rain", "ice"])
def test_substep_bound_covers_golden_columns(golden, name):
    fn, argfn = _KERNELS[name][:2]
    out = fn(**argfn(golden))
    conv = np.asarray(out["converged"])
    assert conv.shape == (golden["st"]["rho"].shape[0],)
    assert conv.all(), (
        f"{name} sedimentation: {(~conv).sum()} of {conv.size} golden "
        "columns did not spend dt within MAX_SEDI_SUBSTEPS")


def test_substep_bound_through_p3_main(golden):
    """Full-pipeline convergence flag (sedimentation runs on the part2
    state inside p3_main) on all 218 golden columns."""
    g, s = golden, golden["s"]
    zcol = g["zcol"]
    out = p3_main(
        g["dt"], True, False, True, False, False,
        g["dry"]["qc"], g["dry"]["nc"], g["dry"]["qr"], g["dry"]["nr"],
        g["dry"]["qi"], g["dry"]["qm"], g["dry"]["ni"], g["dry"]["bm"],
        g["dry"]["qv"], g["th_atm"],
        s["nc_nuceat_tend"], zcol, s["ni_activated"], s["inv_qc_relvar"],
        g["cld_frac_i"], g["cld_frac_l"], g["cld_frac_r"],
        s["p_dry_mid"], g["dz"], s["pseudo_density_dry"], g["inv_exner"],
        calculate_drymmr_from_wetmmr_dp_based(
            s["qv_prev_micro_step"], s["pseudo_density"],
            s["pseudo_density_dry"]),
        s["T_prev_micro_step"], zcol, zcol, zcol, g["tbl"], g["opts"], sed_use_while_loop=False)
    conv = np.asarray(out["sed_converged"])
    assert conv.all(), (
        f"p3_main: {(~conv).sum()} of {conv.size} golden columns did not "
        "spend dt within the sedimentation substep bounds")


# ---------------------------------------------------------------------------
# (i) reverse-mode grad of each kernel: runs, finite, zero NaN
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("name", ["cloud", "rain", "ice"])
def test_sed_kernel_reverse_grad_finite(golden, name):
    scalar, x0 = _scalar_fn(golden, name)
    grad = np.asarray(jax.grad(scalar)(x0))
    n_nan = int(np.isnan(grad).sum())
    n_inf = int(np.isinf(grad).sum())
    assert n_nan == 0, f"{name}: {n_nan}/{grad.size} NaN in reverse grad"
    assert n_inf == 0, f"{name}: {n_inf}/{grad.size} Inf in reverse grad"
    assert np.isfinite(grad).all()
    assert np.linalg.norm(grad) > 0.0  # gradient actually flows


# ---------------------------------------------------------------------------
# (ii) dot-product identity: <grad f, v> == jvp(f; v)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("name", ["cloud", "rain", "ice"])
def test_sed_kernel_grad_jvp_dot_product(golden, name):
    scalar, x0 = _scalar_fn(golden, name)
    v = np.random.default_rng(7).normal(size=x0.shape)
    v /= np.linalg.norm(v)
    grad = np.asarray(jax.grad(scalar)(x0))
    _, tang = jax.jvp(scalar, (x0,), (jnp.asarray(v),))
    dot = float(grad.ravel() @ v.ravel())
    np.testing.assert_allclose(dot, float(tang), rtol=1e-9)


# ---------------------------------------------------------------------------
# (iii) full p3_process_step: reverse grad runs and is NaN-free
# ---------------------------------------------------------------------------
def _process_scalar(golden, wrt, of):
    g, s = golden, golden["s"]
    n = 8

    def scalar(x):
        fields = {k: jnp.asarray(_sub(s[k], n)) for k in
                  ("T_mid", "p_mid", "p_dry_mid", "pseudo_density",
                   "pseudo_density_dry", "cldfrac_tot", "qv", "qc", "nc",
                   "qr", "nr", "qi", "qm", "ni", "bm",
                   "qv_prev_micro_step", "T_prev_micro_step",
                   "nc_nuceat_tend", "ni_activated", "inv_qc_relvar",
                   "precip_liq_surf_mass", "precip_ice_surf_mass")}
        fields[wrt] = x
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
            g["tbl"], g["opts"], sed_use_while_loop=False)
        return sum(jnp.sum(out[k]) for k in of)

    return scalar, jnp.asarray(_sub(s[wrt], n))


def test_p3_process_step_reverse_grad_wrt_qc(golden):
    scalar, x0 = _process_scalar(golden, "qc", ("qc", "qr"))
    grad = np.asarray(jax.grad(scalar)(x0))  # must not raise
    n_nan = int(np.isnan(grad).sum())
    assert n_nan == 0, (
        f"d(sum qc+qr)/d(qc): {n_nan}/{grad.size} NaN in reverse grad")
    assert int(np.isinf(grad).sum()) == 0


def test_p3_process_step_reverse_grad_wrt_qi(golden):
    scalar, x0 = _process_scalar(golden, "qi", ("qi", "qc", "qr"))
    grad = np.asarray(jax.grad(scalar)(x0))  # must not raise
    n_nan = int(np.isnan(grad).sum())
    assert n_nan == 0, (
        f"d(sum qi+qc+qr)/d(qi): {n_nan}/{grad.size} NaN in reverse grad")
    assert int(np.isinf(grad).sum()) == 0
