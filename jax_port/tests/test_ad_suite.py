"""AD health of the ASSEMBLED physics suite (scream_jax.driver ad_mode).

The differentiable path (ScreamPhysics.step(..., ad_mode=True)) must
(1) reproduce the default, bitwise-sacred path exactly on the primal,
(2) yield finite, NaN-free reverse-mode gradients of scalar objectives
through one full suite step w.r.t. T_mid and qc, and (3) stay finite
through a 2-step rollout (output state fed back as input).

State loading follows harness/ad_probe.py: the physics-suite golden
archive provides the fields, the IC file the grid geometry. Skipped if
either is unavailable. Gradient tests run eagerly (no outer jit): the
suite-step XLA compile takes ~20 min while eager reverse mode is ~80 s,
and eager keeps the per-process inner jits cached.
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
DATA = REPO.parent / "e3sm-inputdata" / "atm" / "scream"
IC_FILE = DATA / "init" / "screami_unit_tests_ne2np4L72_20220822.nc"

NCOL = 8
DOY0 = 284 + 45000.0 / 86400.0  # t0 = 2021-10-12-45000 (golden run)


@functools.lru_cache(maxsize=1)
def _load():
    if not GOLDEN.exists():
        pytest.skip(f"golden archive not found: {GOLDEN}")
    if not IC_FILE.exists():
        pytest.skip(f"IC file not found: {IC_FILE}")
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
    phys = ScreamPhysics(
        DATA, geo["hyam"], geo["hybm"], geo["lat"][:NCOL], geo["lon"][:NCOL],
        np.asarray(calculate_dx_from_area(geo["area"][:NCOL],
                                          geo["lat"][:NCOL])),
        mac_mic_subcycles=6, year=2021, spa_col_indices=np.arange(NCOL))
    return phys, s, float(meta["dt"])


def _rollout_objective(phys, s, dt, nsteps, field="T_mid"):
    """sum(field) after an nsteps rollout, as a function of (T0, qc0)."""
    spa_all = [phys.spa_substep_outputs(dt, DOY0 + k * dt / 86400.0,
                                        s["p_mid"])
               for k in range(nsteps)]

    def obj(T0, qc0):
        st = dict(s)
        st["T_mid"] = T0
        st["qc"] = qc0
        for k in range(nsteps):
            st = phys.step(st, dt, k, DOY0 + k * dt / 86400.0,
                           ad_mode=True, spa_outs=spa_all[k])
        return jnp.sum(st[field])
    return obj


def test_ad_mode_primal_matches_default():
    """The traceable path must be a pure refactor of the default path:
    every output field bit-identical (the np round-trips it drops are
    value-preserving, the SPA hoist is exact, and P3's bounded-scan
    sedimentation is bitwise-equal to the while-loop realization)."""
    phys, s, dt = _load()
    out_def = phys.step(dict(s), dt, 0, DOY0)
    out_ad = phys.step(dict(s), dt, 0, DOY0, ad_mode=True)

    assert set(out_ad) == set(out_def)
    bad = {}
    for k in out_def:
        a = np.asarray(out_def[k], dtype=np.float64)
        b = np.asarray(out_ad[k], dtype=np.float64)
        assert a.shape == b.shape, k
        d = np.abs(a - b)
        if d.size and d.max() != 0.0:
            bad[k] = d.max()
    assert not bad, f"ad_mode primal deviates from default path: {bad}"


@pytest.mark.slow
def test_suite_step_reverse_grad_finite():
    """Reverse-mode grad of sum(T_mid out) through ONE full suite step
    (6x [shoc, cld_fraction, spa, p3] + rrtmgp) w.r.t. T_mid and qc:
    must run, be finite, NaN-free, and not identically zero."""
    phys, s, dt = _load()
    obj = _rollout_objective(phys, s, dt, nsteps=1)
    g_T, g_qc = jax.grad(obj, argnums=(0, 1))(
        jnp.asarray(s["T_mid"]), jnp.asarray(s["qc"]))
    for name, g in (("T_mid", g_T), ("qc", g_qc)):
        g = np.asarray(g)
        assert np.isfinite(g).all(), \
            (f"non-finite d/d{name}: NaN={np.isnan(g).sum()} "
             f"Inf={np.isinf(g).sum()} of {g.size}")
        assert np.linalg.norm(g) > 0.0, f"d/d{name} identically zero"


@pytest.mark.slow
def test_rollout2_reverse_grad_finite():
    """Reverse-mode grad through a 2-step rollout (step-1 output state
    fed back as step-2 input; surface/SPA inputs fixed) stays finite."""
    phys, s, dt = _load()
    obj = _rollout_objective(phys, s, dt, nsteps=2)
    g_T, g_qc = jax.grad(obj, argnums=(0, 1))(
        jnp.asarray(s["T_mid"]), jnp.asarray(s["qc"]))
    for name, g in (("T_mid", g_T), ("qc", g_qc)):
        g = np.asarray(g)
        assert np.isfinite(g).all(), \
            (f"non-finite 2-step d/d{name}: NaN={np.isnan(g).sum()} "
             f"Inf={np.isinf(g).sum()} of {g.size}")
        assert np.linalg.norm(g) > 0.0, f"2-step d/d{name} identically zero"
