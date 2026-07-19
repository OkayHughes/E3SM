"""AD health of the SHOC process step on golden states.

Reverse-mode grad of scalar reductions of shoc_process_step must be finite
(no NaN/Inf), must satisfy the dot-product identity against forward-mode
jvp, and jvp must agree with central finite differences. Exercised for two
inputs (T_mid and qv) to cover input-specific code paths.

State loading follows harness/ad_probe.py: the physics-suite golden
archive provides the fields, the IC file the grid geometry. Skipped if
either is unavailable.
"""

import functools
import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import jax  # noqa: E402
import jax.numpy as jnp  # noqa: E402

jax.config.update("jax_enable_x64", True)

REPO = Path(__file__).resolve().parents[2]
GOLDEN = REPO / "jax_port" / "golden" / "physics_suite_218x72_dt1800_2steps.npz"
IC_FILE = (REPO.parent / "e3sm-inputdata" / "atm" / "scream" / "init"
           / "screami_unit_tests_ne2np4L72_20220822.nc")

NCOL = 6
DT = 300.0


@functools.lru_cache(maxsize=1)
def _load():
    if not GOLDEN.exists():
        pytest.skip(f"golden archive not found: {GOLDEN}")
    if not IC_FILE.exists():
        pytest.skip(f"IC file not found: {IC_FILE}")
    netCDF4 = pytest.importorskip("netCDF4")

    z = np.load(GOLDEN)
    meta = json.loads(str(z["__metadata__"]))
    s = {name: np.asarray(z[f"{name}__step0"], dtype=np.float64)[:NCOL]
         for name in meta["fields"]}
    ds = netCDF4.Dataset(IC_FILE)
    geo = {k: np.array(ds[k][:]) for k in ("hyam", "hybm", "lat", "area")}
    ds.close()
    return s, geo


@functools.lru_cache(maxsize=4)
def _scalar_fn(wrt: str):
    """f(x) = sum(T_mid out) of shoc_process_step with x replacing s[wrt].

    Returns (f, x0). Cached so grad/jvp/FD reuse one closure (and its
    compiled trace) per input field.
    """
    s, geo = _load()
    from scream_jax.driver import P0, SHOC_PARAMS
    from scream_jax.foundation.thermo import calculate_dx_from_area
    from scream_jax.shoc import shoc_init
    from scream_jax.shoc.process import shoc_process_step

    npbl = shoc_init(len(geo["hyam"]), 0,
                     P0 * (geo["hyam"] + geo["hybm"]))
    cell_length = jnp.asarray(calculate_dx_from_area(
        geo["area"][:NCOL], geo["lat"][:NCOL]))

    inputs = {
        "T_mid": jnp.asarray(s["T_mid"]),
        "qv": jnp.asarray(s["qv"]),
        "tke": jnp.asarray(s["tke"]),
    }
    assert wrt in inputs

    def f(x):
        kw = dict(inputs)
        kw[wrt] = x
        out = shoc_process_step(
            DT, npbl, cell_length,
            SHOC_PARAMS["lambda_low"], SHOC_PARAMS["lambda_high"],
            SHOC_PARAMS["lambda_slope"], SHOC_PARAMS["lambda_thresh"],
            SHOC_PARAMS["thl2tune"], SHOC_PARAMS["qw2tune"],
            SHOC_PARAMS["qwthl2tune"], SHOC_PARAMS["w2tune"],
            SHOC_PARAMS["length_fac"], SHOC_PARAMS["c_diag_3rd_mom"],
            SHOC_PARAMS["ckh"], SHOC_PARAMS["ckm"], False, False,
            kw["T_mid"], jnp.asarray(s["p_mid"]), jnp.asarray(s["p_int"]),
            jnp.asarray(s["pseudo_density"]), jnp.asarray(s["omega"]),
            jnp.asarray(s["phis"]), jnp.asarray(s["surf_sens_flux"]),
            jnp.asarray(s["surf_evap"]),
            jnp.asarray(s["surf_mom_flux"][:, 0]),
            jnp.asarray(s["surf_mom_flux"][:, 1]),
            kw["qv"], jnp.asarray(s["qc"]), kw["tke"],
            jnp.asarray(s["horiz_winds"][:, 0, :]),
            jnp.asarray(s["horiz_winds"][:, 1, :]),
            jnp.asarray(s["cldfrac_liq"]),
            jnp.asarray(s["sgs_buoy_flux"]),
            jnp.asarray(s["eddy_diff_mom"]))
        return jnp.sum(out["T_mid"])

    return f, np.asarray(inputs[wrt])


@functools.lru_cache(maxsize=4)
def _grad(wrt: str):
    f, x0 = _scalar_fn(wrt)
    return np.asarray(jax.grad(f)(jnp.asarray(x0)))


def _directions(shape, n):
    rng = np.random.default_rng(1234)
    vs = []
    for _ in range(n):
        v = rng.normal(size=shape)
        vs.append(v / np.linalg.norm(v))
    return vs


@pytest.mark.parametrize("wrt", ["T_mid", "qv", "tke"])
def test_grad_finite(wrt):
    """Reverse-mode gradient of sum(T_out) has zero NaN/Inf entries."""
    g = _grad(wrt)
    n_nan = int(np.isnan(g).sum())
    n_inf = int(np.isinf(g).sum())
    assert n_nan == 0, f"grad wrt {wrt}: {n_nan}/{g.size} NaN"
    assert n_inf == 0, f"grad wrt {wrt}: {n_inf}/{g.size} Inf"
    assert np.linalg.norm(g) > 0.0, f"grad wrt {wrt} is identically zero"


@pytest.mark.parametrize("wrt", ["T_mid", "qv", "tke"])
def test_grad_jvp_dot_product_identity(wrt):
    """<grad, v> equals the forward-mode jvp tangent for random v."""
    f, x0 = _scalar_fn(wrt)
    g = _grad(wrt)
    for v in _directions(x0.shape, 3):
        _, tang = jax.jvp(f, (jnp.asarray(x0),), (jnp.asarray(v),))
        tang = float(tang)
        dot = float(g.ravel() @ v.ravel())
        assert tang == pytest.approx(dot, rel=1e-9), \
            f"wrt {wrt}: <grad,v>={dot!r} vs jvp={tang!r}"


# eps (relative to mean |x0|) per input: T_mid at 1e-3; qv at 1e-2 because
# its scale (~2.7e-3) makes eps = 1e-3*scale ~ 2.7e-6, deep in float64
# cancellation noise against an output of ~1.4e5 (FD agreement at 1e-2 is
# 1.3e-4). tke is excluded: the golden state sits exactly on the mintke
# floor, so f is kinked at x0 by check_tke and central FD is meaningless
# there (grad/jvp consistency for tke-adjacent paths is still covered via
# the qv input, which feeds the same tracer/diffusion machinery).
@pytest.mark.parametrize("wrt,eps_rel", [("T_mid", 1e-3), ("qv", 1e-2)])
def test_jvp_vs_central_fd(wrt, eps_rel):
    """Forward-mode jvp matches central finite differences (locks in the
    already-correct forward mode)."""
    f, x0 = _scalar_fn(wrt)
    v = _directions(x0.shape, 1)[0]
    _, tang = jax.jvp(f, (jnp.asarray(x0),), (jnp.asarray(v),))
    tang = float(tang)

    eps = eps_rel * np.abs(x0).mean()
    fp = float(f(jnp.asarray(x0 + eps * v)))
    fm = float(f(jnp.asarray(x0 - eps * v)))
    fd = (fp - fm) / (2.0 * eps)
    assert tang == pytest.approx(fd, rel=1e-3), \
        f"wrt {wrt}: jvp={tang!r} vs central FD={fd!r}"
