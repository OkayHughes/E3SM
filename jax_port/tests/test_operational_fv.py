"""Tests for the operational (HOMME gllfvremap) pg-N layer
(jax_port/pyses_ext/finite_volume_grid_operational.py) against the real
pySEs finite_volume_grid module.

Key checks:
- HOMME's constrained-projection FV->GLL operator (gfr_init_R /
  gfr_f2g_remapd_op, incl. the npi = max(2, nf) intermediate grid) equals
  pySEs' reference fv_to_gll operator -- the claim that lets the
  operational layer reuse the existing reference operators.
- Exact spheremp*dp mass conservation of the dp-weighted remaps and CAAS.
- Shape preservation (bounds) of the mixing-ratio remaps.
- Vector remap: pointwise-tensor consistency and FV->GLL->FV idempotence.
- Full-state drivers: zero physics tendency reproduces the GLL state
  exactly (HOMME's bound-augmentation property); tracer mass bookkeeping.
"""

import os
import sys
from pathlib import Path

import numpy as np
import pytest

PYSES_ROOT = Path(os.environ.get(
    "PYSES_PATH", "/Users/ostensiblyowen/development/python/pyses_07_04_26"))
if not (PYSES_ROOT / "pyses").is_dir():
    pytest.skip("pySEs source tree not available", allow_module_level=True)
sys.path.insert(0, str(PYSES_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pyses_ext"))

from pyses.dynamical_cores.finite_volume_grid import (  # noqa: E402
    _gll_nodes, _interp_low_to_high, _subcell_basis_integrals,
    init_fv_grid, gll_to_fv, reference_fv_operators)

import finite_volume_grid_operational as op  # noqa: E402

NPT, NF, NELEM, NLEV = 4, 2, 6, 10


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------
def make_grids(nf=NF, seed=0):
    rng = np.random.default_rng(seed)
    metdet = 1.0 + 0.3 * rng.uniform(-1, 1, (NELEM, NPT, NPT))
    c2p = (np.eye(2)[None, None, None]
           + 0.4 * rng.uniform(-1, 1, (NELEM, NPT, NPT, 2, 2)))
    p2c = np.linalg.inv(c2p)
    h_grid = {
        "metric_determinant": metdet,
        "physical_coords": rng.uniform(-1, 1, (NELEM, NPT, NPT, 2)),
        "contra_to_physical": c2p,
        "physical_to_contra": p2c,
    }
    dims = {"npt": NPT, "num_elem": NELEM}
    fv_grid = init_fv_grid(h_grid, dims, nf=nf)
    fv_grid = op.extend_fv_grid_operational(fv_grid, h_grid)

    ai = np.linspace(0.02, 0.0, NLEV + 1) ** 1.5 + 1e-3
    bi = np.linspace(0.0, 1.0, NLEV + 1) ** 1.2
    v_grid = {
        "hybrid_a_i": ai, "hybrid_b_i": bi,
        "hybrid_a_m": 0.5 * (ai[1:] + ai[:-1]),
        "hybrid_b_m": 0.5 * (bi[1:] + bi[:-1]),
        "reference_surface_mass": 1.0e5,
    }
    return h_grid, fv_grid, v_grid, rng


def elem_mass(q, dp, spheremp):
    """Per-element spheremp*dp-weighted tracer mass."""
    arr = np.asarray(q) * np.asarray(dp)
    w = np.asarray(spheremp)
    if arr.ndim == 4:
        w = w[..., None]
    return (arr * w).reshape(arr.shape[0], -1).sum(1)


# ---------------------------------------------------------------------------
# operator equivalence: HOMME constrained projection == pySEs reference
# ---------------------------------------------------------------------------
def homme_f2g_1d(npt, nf):
    """1-D FV->GLL operator built the HOMME way (gfr_init_R /
    gfr_f2g_remapd_op): mass-constrained L2 projection onto the
    npi = max(2, nf) GLL basis, then interpolation to the np basis.
    The 2-D operator is this matrix tensor-produced with itself."""
    npi = max(2, nf)
    gll_i = _gll_nodes(npi)
    w_i = op._gll_weights(npi)
    m1, _ = _subcell_basis_integrals(gll_i, nf)         # (nf, npi)
    A = m1.T                                            # M_sgf, (npi, nf)
    wff1 = 2.0 / nf                                     # M_ff = wff1 * I
    S = A.T @ (A / w_i[:, None])
    g_npi = (A @ np.linalg.solve(S, wff1 * np.eye(nf))) / w_i[:, None]
    interp = _interp_low_to_high(gll_i, _gll_nodes(npt))  # (npt, npi)
    return interp @ g_npi                               # (npt, nf)


@pytest.mark.parametrize("nf", [2, 3, 4])
def test_constrained_projection_equals_reference(nf):
    ref = reference_fv_operators(NPT, nf)["fv_to_gll"]
    homme = homme_f2g_1d(NPT, nf)
    np.testing.assert_allclose(homme, ref, atol=1e-12)


def test_constrained_projection_pg1():
    # pySEs' reference operators don't tabulate a 1-node basis; HOMME's
    # npi=2 least-squares branch reduces analytically to the constant
    # (piecewise-constant FV -> constant GLL) operator for pg1.
    homme = homme_f2g_1d(NPT, 1)
    np.testing.assert_allclose(homme, np.ones((NPT, 1)), atol=1e-13)


def test_gll_weights():
    # against pySEs' tabulated 4-point weights
    np.testing.assert_allclose(op._gll_weights(4),
                               [1 / 6, 5 / 6, 5 / 6, 1 / 6], atol=1e-14)
    np.testing.assert_allclose(op._gll_weights(2), [1.0, 1.0], atol=1e-14)


# ---------------------------------------------------------------------------
# dp-weighted remaps
# ---------------------------------------------------------------------------
def test_dp_weighted_mass_conservation():
    h_grid, fv_grid, v_grid, rng = make_grids()
    ps = 1e5 * (1 + 0.05 * rng.uniform(-1, 1, (NELEM, NPT, NPT)))
    dp_g = op.surface_mass_to_d_mass(ps, v_grid)
    dp_f = op.calc_dp_fv(gll_to_fv(ps, fv_grid), v_grid)
    q = rng.uniform(0.0, 1e-2, (NELEM, NPT, NPT, NLEV))

    q_f = op.gll_to_fv_dp(q, dp_g, dp_f, fv_grid)
    m_g = elem_mass(q, dp_g, fv_grid["spheremp"])
    m_f = elem_mass(q_f, dp_f, fv_grid["spheremp_fv"])
    np.testing.assert_allclose(m_f, m_g, rtol=1e-13)

    q_b = op.fv_to_gll_dp(q_f, dp_f, dp_g, fv_grid)
    m_b = elem_mass(q_b, dp_g, fv_grid["spheremp"])
    np.testing.assert_allclose(m_b, m_f, rtol=1e-13)


def test_dp_weighted_constant_is_exact():
    h_grid, fv_grid, v_grid, rng = make_grids()
    ps = 1e5 * (1 + 0.05 * rng.uniform(-1, 1, (NELEM, NPT, NPT)))
    dp_g = op.surface_mass_to_d_mass(ps, v_grid)
    dp_f = op.calc_dp_fv(gll_to_fv(ps, fv_grid), v_grid)
    q = np.full((NELEM, NPT, NPT, NLEV), 3.25e-3)
    q_f = op.gll_to_fv_dp(q, dp_g, dp_f, fv_grid)
    # constant q: remap(dp*q)/dp_f = q * remap(dp)/dp_f; dp_f is the
    # hydrostatic reconstruction of remapped ps, which by linearity of the
    # remap and the hybrid formula equals remap(dp) exactly
    np.testing.assert_allclose(q_f, 3.25e-3, rtol=1e-12)


# ---------------------------------------------------------------------------
# CAAS limiter
# ---------------------------------------------------------------------------
def test_caas_noop_within_bounds():
    rng = np.random.default_rng(1)
    x = rng.uniform(0.2, 0.8, (3, 9, 4))
    c = rng.uniform(0.5, 2.0, (3, 9, 4))
    out = op._caas_clip_and_sum(x, c, np.zeros((3, 4)), np.ones((3, 4)))
    np.testing.assert_array_equal(out, x)


def test_caas_clips_and_conserves():
    rng = np.random.default_rng(2)
    x = rng.uniform(0.0, 1.0, (5, 16))
    x[:, 0] = 2.5    # violates qmax
    x[:, 1] = -0.5   # violates qmin
    c = rng.uniform(0.5, 2.0, (5, 16))
    qmin, qmax = np.full(5, 0.0), np.full(5, 1.0)
    out = op._caas_clip_and_sum(x, c, qmin, qmax)
    assert (out >= -1e-15).all() and (out <= 1 + 1e-15).all()
    np.testing.assert_allclose((c * out).sum(1), (c * x).sum(1), rtol=1e-13)


def test_caas_infeasible_conserves_mass():
    # mean above qmax: bound must be relaxed, mass conserved (HOMME choice)
    x = np.full((1, 4), 3.0)
    c = np.ones((1, 4))
    out = op._caas_clip_and_sum(x, c, np.zeros(1), np.ones(1))
    np.testing.assert_allclose((c * out).sum(1), 12.0, rtol=1e-14)


def test_mixing_ratio_remap_bounds_and_mass():
    h_grid, fv_grid, v_grid, rng = make_grids()
    ps = 1e5 * (1 + 0.05 * rng.uniform(-1, 1, (NELEM, NPT, NPT)))
    dp_g = op.surface_mass_to_d_mass(ps, v_grid)
    dp_f = op.calc_dp_fv(gll_to_fv(ps, fv_grid), v_grid)
    q = rng.uniform(0.0, 1e-2, (NELEM, NPT, NPT, NLEV))
    q[:, 0, 0, :] = 5e-2   # sharp in-element spike

    q_f = op.gll_to_fv_mixing_ratio(q, dp_g, dp_f, fv_grid)
    qmin = q.reshape(NELEM, -1, NLEV).min(1)
    qmax = q.reshape(NELEM, -1, NLEV).max(1)
    assert (q_f.reshape(NELEM, -1, NLEV) >= qmin[:, None] - 1e-15).all()
    assert (q_f.reshape(NELEM, -1, NLEV) <= qmax[:, None] + 1e-15).all()
    np.testing.assert_allclose(
        elem_mass(q_f, dp_f, fv_grid["spheremp_fv"]),
        elem_mass(q, dp_g, fv_grid["spheremp"]), rtol=1e-13)


# ---------------------------------------------------------------------------
# vector remap
# ---------------------------------------------------------------------------
def test_vector_remap_identity_tensors_reduce_to_scalar():
    h_grid, fv_grid, v_grid, rng = make_grids()
    eye = np.broadcast_to(np.eye(2), (NELEM, NPT, NPT, 2, 2)).copy()
    h_id = dict(h_grid, contra_to_physical=eye, physical_to_contra=eye)
    fv_id = init_fv_grid(h_id, {"npt": NPT, "num_elem": NELEM}, nf=NF)
    fv_id = op.extend_fv_grid_operational(fv_id, h_id)
    u = rng.uniform(-30, 30, (NELEM, NPT, NPT, NLEV, 2))
    got = op.gll_to_fv_vector(u, h_id, fv_id)
    want = np.stack([np.asarray(op._ref_g2f(u[..., d], fv_id))
                     for d in range(2)], axis=-1)
    np.testing.assert_allclose(got, want, atol=1e-12)


def test_vector_remap_fv_gll_fv_idempotent():
    h_grid, fv_grid, v_grid, rng = make_grids()
    u_fv = rng.uniform(-30, 30, (NELEM, NF, NF, NLEV, 2))
    u_gll = op.fv_to_gll_vector(u_fv, h_grid, fv_grid)
    back = op.gll_to_fv_vector(u_gll, h_grid, fv_grid)
    np.testing.assert_allclose(back, u_fv, rtol=1e-11, atol=1e-11)


# ---------------------------------------------------------------------------
# full-state drivers
# ---------------------------------------------------------------------------
def make_state(rng, ps_uniform=False):
    if ps_uniform:
        ps = np.full((NELEM, NPT, NPT), 1.0e5)
    else:
        ps = 1e5 * (1 + 0.05 * rng.uniform(-1, 1, (NELEM, NPT, NPT)))
    T = 250.0 + 40.0 * rng.uniform(-1, 1, (NELEM, NPT, NPT, NLEV))
    uv = rng.uniform(-30, 30, (NELEM, NPT, NPT, NLEV, 2))
    omega = rng.uniform(-1, 1, (NELEM, NPT, NPT, NLEV))
    phis = rng.uniform(0, 3e4, (NELEM, NPT, NPT))
    q = {"qv": rng.uniform(1e-6, 2e-2, (NELEM, NPT, NPT, NLEV)),
         "qc": rng.uniform(0.0, 1e-4, (NELEM, NPT, NPT, NLEV))}
    return ps, phis, T, uv, omega, q


def test_dyn_to_fv_phys_uniform_T_exact():
    h_grid, fv_grid, v_grid, rng = make_grids()
    ps, phis, T, uv, omega, q = make_state(rng, ps_uniform=True)
    T[:] = 287.5
    dp_g = op.surface_mass_to_d_mass(ps, v_grid)
    out = op.dyn_to_fv_phys(ps, phis, T, uv, omega, q, dp_g,
                            h_grid, fv_grid, v_grid)
    # uniform ps => uniform p on both grids => theta-form roundtrip exact
    np.testing.assert_allclose(out["T"], 287.5, rtol=1e-12)
    np.testing.assert_allclose(out["ps"], 1.0e5, rtol=1e-13)
    np.testing.assert_allclose(np.asarray(out["dp"]),
                               np.asarray(dp_g[:, :NF, :NF, :]), rtol=1e-12)


def test_dyn_to_fv_phys_phis_limited():
    h_grid, fv_grid, v_grid, rng = make_grids()
    ps, phis, T, uv, omega, q = make_state(rng)
    dp_g = op.surface_mass_to_d_mass(ps, v_grid)
    out = op.dyn_to_fv_phys(ps, phis, T, uv, omega, q, dp_g,
                            h_grid, fv_grid, v_grid)
    lo = phis.reshape(NELEM, -1).min(1)[:, None, None]
    hi = phis.reshape(NELEM, -1).max(1)[:, None, None]
    assert (np.asarray(out["phis"]) >= lo - 1e-9).all()
    assert (np.asarray(out["phis"]) <= hi + 1e-9).all()
    for name in q:
        qf = np.asarray(out["q"][name]).reshape(NELEM, -1, NLEV)
        assert (qf >= q[name].reshape(NELEM, -1, NLEV).min(1)[:, None]
                - 1e-15).all()


def test_fv_phys_to_dyn_zero_tendency_is_identity():
    h_grid, fv_grid, v_grid, rng = make_grids()
    ps, phis, T, uv, omega, q = make_state(rng)
    dp_g = op.surface_mass_to_d_mass(ps, v_grid)
    fwd = op.dyn_to_fv_phys(ps, phis, T, uv, omega, q, dp_g,
                            h_grid, fv_grid, v_grid)
    out = op.fv_phys_to_dyn(np.zeros_like(np.asarray(fwd["T"])),
                            np.zeros_like(np.asarray(fwd["uv"])),
                            fwd["q"], ps, q, dp_g, h_grid, fv_grid, v_grid)
    np.testing.assert_array_equal(np.asarray(out["FT"]), 0.0)
    np.testing.assert_array_equal(np.asarray(out["FM"]), 0.0)
    for name in q:
        np.testing.assert_array_equal(np.asarray(out["FQ"][name]), q[name])


def test_fv_phys_to_dyn_tracer_mass_bookkeeping():
    h_grid, fv_grid, v_grid, rng = make_grids()
    ps, phis, T, uv, omega, q = make_state(rng)
    dp_g = op.surface_mass_to_d_mass(ps, v_grid)
    dp_f = op.calc_dp_fv(gll_to_fv(ps, fv_grid), v_grid)
    fwd = op.dyn_to_fv_phys(ps, phis, T, uv, omega, q, dp_g,
                            h_grid, fv_grid, v_grid)
    # physics adds a bounded increment on the FV grid
    q1_fv = {name: np.asarray(fwd["q"][name])
             * (1 + 0.1 * rng.uniform(-1, 1, (NELEM, NF, NF, NLEV)))
             for name in q}
    out = op.fv_phys_to_dyn(np.zeros((NELEM, NF, NF, NLEV)),
                            np.zeros((NELEM, NF, NF, NLEV, 2)),
                            q1_fv, ps, q, dp_g, h_grid, fv_grid, v_grid)
    for name in q:
        # element mass of FQ == GLL mass + FV-grid physics mass increment
        want = (elem_mass(q[name], dp_g, fv_grid["spheremp"])
                + elem_mass(q1_fv[name] - np.asarray(fwd["q"][name]),
                            dp_f, fv_grid["spheremp_fv"]))
        got = elem_mass(np.asarray(out["FQ"][name]), dp_g,
                        fv_grid["spheremp"])
        np.testing.assert_allclose(got, want, rtol=1e-12)
        # bounds: within [min(FV state, GLL state), max(...)] per element
        lo = np.minimum(
            q1_fv[name].reshape(NELEM, -1, NLEV).min(1),
            q[name].reshape(NELEM, -1, NLEV).min(1))[:, None]
        hi = np.maximum(
            q1_fv[name].reshape(NELEM, -1, NLEV).max(1),
            q[name].reshape(NELEM, -1, NLEV).max(1))[:, None]
        fq = np.asarray(out["FQ"][name]).reshape(NELEM, -1, NLEV)
        assert (fq >= lo - 1e-14).all() and (fq <= hi + 1e-14).all()


def test_neighbor_minmax_hook_is_used():
    h_grid, fv_grid, v_grid, rng = make_grids()
    ps, phis, T, uv, omega, q = make_state(rng)
    dp_g = op.surface_mass_to_d_mass(ps, v_grid)
    fwd = op.dyn_to_fv_phys(ps, phis, T, uv, omega, q, dp_g,
                            h_grid, fv_grid, v_grid)
    calls = []

    def widen(qmin, qmax):
        calls.append(1)
        return qmin - 1.0, qmax + 1.0

    op.fv_phys_to_dyn(np.zeros((NELEM, NF, NF, NLEV)),
                      np.zeros((NELEM, NF, NF, NLEV, 2)),
                      fwd["q"], ps, q, dp_g, h_grid, fv_grid, v_grid,
                      neighbor_minmax=widen)
    assert len(calls) == len(q)
