"""AD hardening tests for the RRTMGP process step.

The radiation step must be traceable under JAX transforms: jax.jvp and
jax.grad of scalar reductions of its outputs must run on golden states,
produce finite derivatives (no NaN/Inf), satisfy the forward/reverse
dot-product identity, and match finite differences on the smooth
d(sum LW_flux_up)/d(T_mid) path.

Uses the same golden-state loading pattern as harness/ad_probe.py and
skips when the golden archive or input data are absent. ncol is kept
small for runtime; the expensive jvp/grad evaluations are shared via a
module-scoped fixture.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import jax  # noqa: E402
import jax.numpy as jnp  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
GOLDEN = REPO / "jax_port" / "golden" / "physics_suite_218x72_dt1800_2steps.npz"
DATA = REPO.parent / "e3sm-inputdata" / "atm" / "scream"
IC_FILE = DATA / "init" / "screami_unit_tests_ne2np4L72_20220822.nc"

NCOL = 6


@pytest.fixture(scope="module")
def ad_setup():
    """Golden state + a scalar-valued radiation step f(T_mid) =
    sum(LW_flux_up), plus its grad and jvp along fixed directions."""
    if not GOLDEN.exists() or not IC_FILE.exists():
        pytest.skip("golden archive or input data not available")
    nc4 = pytest.importorskip("netCDF4")

    z = np.load(GOLDEN)
    meta = json.loads(str(z["__metadata__"]))
    s = {name: np.asarray(z[f"{name}__step0"], dtype=np.float64)[:NCOL]
         for name in meta["fields"]}
    ds = nc4.Dataset(IC_FILE)
    geo = {k: np.array(ds[k][:]) for k in ("hyam", "hybm", "lat", "lon",
                                           "area")}
    ds.close()

    from scream_jax.driver import ScreamPhysics
    from scream_jax.foundation.thermo import calculate_dx_from_area
    from scream_jax.rrtmgp.process import rrtmgp_process_step

    phys = ScreamPhysics(DATA, geo["hyam"], geo["hybm"], geo["lat"][:NCOL],
                         geo["lon"][:NCOL],
                         np.asarray(calculate_dx_from_area(
                             geo["area"][:NCOL], geo["lat"][:NCOL])))

    def rad_lw_up(T):
        out = rrtmgp_process_step(
            phys.kd_sw, phys.kd_lw, phys.co_sw, phys.co_lw,
            phys.rrtmgp_params, 1800.0, 0, 2021, 285.5,
            jnp.asarray(phys.lat), jnp.asarray(phys.lon),
            T, jnp.asarray(s["p_mid"]), jnp.asarray(s["p_int"]),
            jnp.asarray(s["pseudo_density"]),
            jnp.asarray(s["sfc_alb_dir_vis"]),
            jnp.asarray(s["sfc_alb_dir_nir"]),
            jnp.asarray(s["sfc_alb_dif_vis"]),
            jnp.asarray(s["sfc_alb_dif_nir"]),
            jnp.asarray(s["qv"]), jnp.asarray(s["qc"]),
            jnp.asarray(s["nc"]), jnp.asarray(s["qi"]),
            jnp.asarray(s["cldfrac_tot"]),
            jnp.asarray(s["eff_radius_qc"]),
            jnp.asarray(s["eff_radius_qi"]),
            jnp.asarray(s["surf_lw_flux_up"]),
            jnp.asarray(s["o3_volume_mix_ratio"]),
            jnp.asarray(s["rad_heating_pdel"]),
            aero_tau_sw=jnp.asarray(s["aero_tau_sw"]),
            aero_ssa_sw=jnp.asarray(s["aero_ssa_sw"]),
            aero_g_sw=jnp.asarray(s["aero_g_sw"]),
            aero_tau_lw=jnp.asarray(s["aero_tau_lw"]))
        return out["LW_flux_up"]

    def scalar(T):
        return jnp.sum(rad_lw_up(T))

    T0 = jnp.asarray(s["T_mid"])
    rng = np.random.default_rng(0)
    dirs = []
    for _ in range(2):
        v = rng.normal(size=T0.shape)
        dirs.append(v / np.linalg.norm(v))

    primal0 = float(scalar(T0))
    jvps = []
    for v in dirs:
        _, tang = jax.jvp(scalar, (T0,), (jnp.asarray(v),))
        jvps.append(float(tang))
    grad = np.asarray(jax.grad(scalar)(T0))

    return {"scalar": scalar, "T0": T0, "dirs": dirs,
            "primal0": primal0, "jvps": jvps, "grad": grad}


def test_jvp_runs_and_is_finite(ad_setup):
    assert np.isfinite(ad_setup["primal0"])
    assert ad_setup["primal0"] > 0.0  # LW up flux sum is positive
    for tang in ad_setup["jvps"]:
        assert np.isfinite(tang)


def test_grad_runs_no_nan_inf(ad_setup):
    g = ad_setup["grad"]
    assert g.shape == ad_setup["T0"].shape
    n_nan = int(np.isnan(g).sum())
    n_inf = int(np.isinf(g).sum())
    assert n_nan == 0, f"gradient has {n_nan}/{g.size} NaNs"
    assert n_inf == 0, f"gradient has {n_inf}/{g.size} Infs"
    # LW up should actually respond to temperature somewhere
    assert np.linalg.norm(g) > 0.0


def test_grad_jvp_dot_product_identity(ad_setup):
    """<grad f, v> == jvp(f; v) for random directions."""
    g = ad_setup["grad"]
    for v, tang in zip(ad_setup["dirs"], ad_setup["jvps"]):
        dot = float(g.ravel() @ np.asarray(v).ravel())
        assert dot == pytest.approx(tang, rel=1e-10), \
            f"<grad,v>={dot!r} vs jvp={tang!r}"


def test_jvp_matches_finite_difference(ad_setup):
    """Central FD of d(sum LW_flux_up)/d(T_mid) at eps=1e-3 relative;
    LW is smooth in T away from branch points."""
    scalar = ad_setup["scalar"]
    T0 = ad_setup["T0"]
    v = ad_setup["dirs"][0]
    tang = ad_setup["jvps"][0]

    eps = 1e-3 * float(jnp.abs(T0).mean())
    fp = float(scalar(T0 + eps * jnp.asarray(v)))
    fm = float(scalar(T0 - eps * jnp.asarray(v)))
    fd = (fp - fm) / (2.0 * eps)
    rel = abs(fd - tang) / max(abs(tang), 1e-30)
    assert rel < 2e-3, f"fd={fd!r} jvp={tang!r} rel={rel:.3e}"
