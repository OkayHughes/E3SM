"""Operational (HOMME ``gllfvremap``) extension of the pg-N physics grid.

This module adds the *application layer* that SCREAM/E3SM uses operationally
on top of the element-local remap operators already implemented in
``pyses/dynamical_cores/finite_volume_grid.py``.  It is written to be appended
to that file: every function below only uses symbols defined there (plus
``pyses.dynamical_cores.mass_coordinate``), 2-space indentation and the same
array conventions.  When merging, delete the import block at the top.

Source of truth: ``components/homme/src/share/gllfvremap_mod.F90`` (E3SM,
A.M. Bradley 2019-2020), the remap used by EAMxx/SCREAM for pg2 physics
(``eamxx_homme_fv_phys.cpp`` -> Hommexx ``GllFvRemap``, whose operators are
built by that Fortran module).  The *reference operators* there coincide with
:func:`reference_fv_operators` for ``nf >= 2``: HOMME builds FV -> GLL as the
mass-constrained L2 projection onto an intermediate ``npi = max(2, nf)`` GLL
basis (``gfr_init_R``/``gfr_f2g_remapd_op``), which for ``nf == npi`` reduces
exactly to the paper operator ``I (A^{p->f} I)^{-1}`` used here (see the
regression test ``test_constrained_projection_equals_reference``).  What does
NOT coincide -- and what this module adds -- is how the operators are applied
to a model state:

* **dp-weighted remap of mixing ratios** (``gfr_g2f_scalar_dp`` /
  ``gfr_f2g_scalar_dp``): remap ``dp*q`` as a density, divide by the target
  grid's ``dp``;
* **shape preservation** (``gfr_g2f_mixing_ratio``, ``limiter_clip_and_sum``):
  the CAAS limiter (Alg 3.1 of doi:10.1137/18M1165414) restores per-element
  extrema bounds after remap while conserving the ``spheremp*dp``-weighted
  mass, applied to the FULL tracer state (not the tendency);
* **theta-form temperature remap**: T is converted to ``T*(p0/p)^kappa``,
  remapped dp-weighted, and converted back with the target grid's own
  (hydrostatically reconstructed) pressure;
* **hydrostatic FV pressure**: ``dp_fv`` is rebuilt from remapped ``ps`` via
  the hybrid coefficients (``calc_dp_fv``) rather than remapped directly --
  "this loop is essentially how CAM computes pdel, so we must use it, too";
* **vector remap through the element map** (``gfr_g2f_vector`` /
  ``gfr_f2g_vector``): winds are converted physical -> contravariant
  (reference-element) components, remapped *unweighted*, and converted back
  with tensors evaluated on the target grid, avoiding steep component
  gradients at the poles;
* **tendency/state asymmetry** (``gfr_fv_phys_to_dyn``): T and winds return
  to the dynamics grid as *tendencies*; tracers return as *full states*,
  limited against halo-exchanged FV bounds augmented with the pre-physics GLL
  bounds (so a zero physics tendency reproduces the GLL state exactly).

Deliberate divergences from HOMME, all element-local by construction:

* The FV-point map tensors ``contra_to_physical_fv`` are obtained by
  remapping the GLL tensor entries and rescaling their determinant to the
  remapped metric determinant, instead of evaluating the analytic
  element->sphere map at the FV cell centres (``gfr_init_Dmap``); this is a
  within-element interpolation of a smooth tensor and keeps the module
  independent of pySEs' mesh-generation internals.  Evaluate the map exactly
  there if/when the FV grid init gains access to it.
* Halo exchange of the limiter bounds (``gfr_f2g_mixing_ratios_he``) is a
  caller-supplied hook (``neighbor_minmax``); the default is element-local
  bounds.  A conforming hook can be built from
  ``pyses.operations_2d.local_assembly.minmax_scalar``.
* The pg1 order-of-accuracy boost (``gfr_pg1_init``) is not implemented
  (``nf == 1`` still works, without the boost).
* No DSS is performed here (HOMME's ``gfr_f2g_dss`` runs afterwards, too):
  apply pySEs' assembly to the returned FT/FM/FQ to make them continuous.
"""
import numpy as np

# --- standalone-import shim -------------------------------------------------
# DELETE this block when merging into pyses/dynamical_cores/finite_volume_grid
# (every name below is already in that module's scope; ``surface_mass_to_d_mass``
# then needs `from .mass_coordinate import surface_mass_to_d_mass` at the top).
from pyses.dynamical_cores.finite_volume_grid import (  # noqa: F401
    bnp, device_wrapper, unwrap, _gll_nodes, _lagrange_legendre_coeffs,
    _subcell_basis_integrals, _safe_divide, reference_fv_operators,
    init_fv_grid, gll_to_fv, fv_to_gll, fv_field_to_columns,
    columns_to_fv_field)
from pyses.dynamical_cores.mass_coordinate import surface_mass_to_d_mass
# --- end shim ----------------------------------------------------------------

#: E3SM/EAMxx dry-air kappa = R_dry / cp_dry.
KAPPA_DRY = 287.042 / 1004.64


def _gll_weights(n):
  """1-D GLL quadrature weights for the :func:`_gll_nodes` ordering.

  ``w_i = 2 / (n (n-1) P_{n-1}(x_i)^2)`` -- exact for any ``n >= 2``, so the
  2-point rule (needed for the panel basis) works without a table entry.
  """
  x = _gll_nodes(n)
  c = np.zeros(n)
  c[n - 1] = 1.0
  pnm1 = np.polynomial.legendre.legval(x, c)
  return 2.0 / (n * (n - 1) * pnm1 ** 2)


# ---------------------------------------------------------------------------
# Operational FV grid struct
# ---------------------------------------------------------------------------
def extend_fv_grid_operational(fv_grid, h_grid):
  """Add the operational (HOMME ``gllfvremap``) data to an FV grid struct.

  Parameters
  ----------
  fv_grid : dict
      Output of :func:`init_fv_grid`.
  h_grid : SpectralElementGrid
      Must additionally contain ``"contra_to_physical"`` and
      ``"physical_to_contra"`` ``(E, np, np, 2, 2)`` (the vector remap
      tensors; see :func:`gll_to_fv_vector`).

  Returns
  -------
  fv_grid : dict
      The same struct, extended with
      ``"gll_weights_2d"`` ``(np, np)`` (GLL quadrature weights, HOMME
      ``w_gg``), ``"w_ff"`` (scalar FV reference cell area ``(2/nf)^2``,
      HOMME ``w_ff``), ``"spheremp"`` ``(E, np, np)`` and ``"spheremp_fv"``
      ``(E, nf, nf)`` (the metric-weighted quadrature weights HOMME uses as
      limiter masses), and ``"contra_to_physical_fv"`` /
      ``"physical_to_contra_fv"`` ``(E, nf, nf, 2, 2)``.
  """
  npt, nf = fv_grid["npt"], fv_grid["nf"]
  w1 = _gll_weights(npt)
  w_gg = np.outer(w1, w1)                                        # (np, np)
  w_ff = (2.0 / nf) ** 2

  metdet = np.asarray(unwrap(fv_grid["metric_determinant"]))     # (E, np, np)
  metdet_fv = np.asarray(unwrap(fv_grid["metric_determinant_fv"]))

  # FV-point map tensors: remap the (smooth, element-local) GLL tensor
  # entries, then rescale so det(D_fv) equals the remapped determinant --
  # HOMME's normalization det(D) == metdet (gfr_init_Dmap), expressed in
  # whatever normalization h_grid's tensors use.
  c2p = np.asarray(unwrap(h_grid["contra_to_physical"]))         # (E,np,np,2,2)
  ops = reference_fv_operators(npt, nf)
  a2f = ops["gll_to_fv"]                                         # (nf, np)
  c2p_fv = np.einsum("ci,dj,eijab->ecdab", a2f, a2f, c2p)
  det_gll = c2p[..., 0, 0] * c2p[..., 1, 1] - c2p[..., 0, 1] * c2p[..., 1, 0]
  det_tgt = np.einsum("ci,dj,eij,eij->ecd", a2f, a2f, metdet, det_gll)
  det_tgt = _safe_divide(det_tgt, np.asarray(metdet_fv))
  det_fv = (c2p_fv[..., 0, 0] * c2p_fv[..., 1, 1]
            - c2p_fv[..., 0, 1] * c2p_fv[..., 1, 0])
  scale = np.sqrt(np.abs(_safe_divide(det_tgt, det_fv)))
  c2p_fv = c2p_fv * scale[..., None, None]
  det_fv = det_fv * scale ** 2
  inv_det = np.asarray(_safe_divide(np.ones_like(det_fv), det_fv))
  p2c_fv = np.empty_like(c2p_fv)                                 # (g, s) layout
  p2c_fv[..., 0, 0] = c2p_fv[..., 1, 1] * inv_det
  p2c_fv[..., 0, 1] = -c2p_fv[..., 0, 1] * inv_det
  p2c_fv[..., 1, 0] = -c2p_fv[..., 1, 0] * inv_det
  p2c_fv[..., 1, 1] = c2p_fv[..., 0, 0] * inv_det

  fv_grid["gll_weights_2d"] = device_wrapper(w_gg)
  fv_grid["w_ff"] = float(w_ff)
  fv_grid["spheremp"] = device_wrapper(w_gg[None] * metdet,
                                       elem_sharding_axis=0)
  fv_grid["spheremp_fv"] = device_wrapper(w_ff * metdet_fv,
                                          elem_sharding_axis=0)
  fv_grid["contra_to_physical_fv"] = device_wrapper(c2p_fv,
                                                    elem_sharding_axis=0)
  fv_grid["physical_to_contra_fv"] = device_wrapper(p2c_fv,
                                                    elem_sharding_axis=0)
  return fv_grid


# ---------------------------------------------------------------------------
# CAAS limiter (HOMME limiter1_clip_and_sum)
# ---------------------------------------------------------------------------
def _caas_clip_and_sum(x, c, qmin, qmax):
  """Clip-and-sum limiter (CAAS, Alg 3.1 of doi:10.1137/18M1165414).

  Solves ``min_x* ||c (x - x*)||_1  st  c'x* = c'x,  qmin <= x* <= qmax``
  over the cell axis (axis 1).  ``x``/``c`` are ``(E, n[, K])``;
  ``qmin``/``qmax`` are ``(E[, K])``.  When the problem is infeasible the
  violated bound is relaxed to the mean, conserving mass (as in HOMME).
  """
  sumc = bnp.sum(c, axis=1)
  mass = bnp.sum(c * x, axis=1)
  mean = _safe_divide(mass, sumc)
  qmin = bnp.where(mass < qmin * sumc, mean, qmin)
  qmax = bnp.where(mass > qmax * sumc, mean, qmax)
  lo = bnp.expand_dims(qmin, 1)
  hi = bnp.expand_dims(qmax, 1)
  xc = bnp.clip(x, lo, hi)
  addmass = mass - bnp.sum(c * xc, axis=1)
  v = bnp.where(bnp.expand_dims(addmass, 1) > 0.0, hi - xc, xc - lo)
  den = bnp.sum(v * c, axis=1)
  fac = bnp.where(den > 0.0, _safe_divide(addmass, den), 0.0)
  return xc + bnp.expand_dims(fac, 1) * v


def _flatten_cells(field):
  """``(E, a, a, ...) -> (E, a*a, ...)`` (in-element cell axis for CAAS)."""
  return bnp.reshape(field,
                     (field.shape[0], -1) + tuple(field.shape[3:]))


# ---------------------------------------------------------------------------
# dp-weighted scalar and mixing-ratio remaps
# ---------------------------------------------------------------------------
def gll_to_fv_dp(q, dp_gll, dp_fv, fv_grid):
  """dp-weighted GLL -> FV remap of a mixing ratio (``gfr_g2f_scalar_dp``).

  ``dp*q`` is remapped as a density and divided by the FV ``dp``; conserves
  ``sum(spheremp*dp*q)`` exactly for any (positive) ``dp`` pair.  Shapes:
  ``q``/``dp_gll`` ``(E, np, np[, K])``, ``dp_fv`` ``(E, nf, nf[, K])``.
  """
  return _safe_divide(gll_to_fv(q * dp_gll, fv_grid), dp_fv)


def fv_to_gll_dp(q_fv, dp_fv, dp_gll, fv_grid):
  """dp-weighted FV -> GLL remap of a mixing ratio (``gfr_f2g_scalar_dp``)."""
  return _safe_divide(fv_to_gll(q_fv * dp_fv, fv_grid), dp_gll)


def gll_to_fv_mixing_ratio(q, dp_gll, dp_fv, fv_grid):
  """Shape-preserving dp-weighted GLL -> FV remap (``gfr_g2f_mixing_ratio``).

  After the dp-weighted remap, CAAS restores the per-element(-level) GLL
  extrema bounds while conserving the ``spheremp_fv*dp_fv``-weighted mass.
  """
  q_fv = gll_to_fv_dp(q, dp_gll, dp_fv, fv_grid)
  qmin = bnp.min(_flatten_cells(q), axis=1)
  qmax = bnp.max(_flatten_cells(q), axis=1)
  c = _flatten_cells(dp_fv * _with_levels(fv_grid["spheremp_fv"], dp_fv))
  out = _caas_clip_and_sum(_flatten_cells(q_fv), c, qmin, qmax)
  return bnp.reshape(out, q_fv.shape)


def gll_to_fv_limited(g, fv_grid):
  """Limited (but not dp-weighted) GLL -> FV scalar remap
  (``gfr_g2f_scalar_and_limit``) -- HOMME uses this for ``phis``."""
  f = gll_to_fv(g, fv_grid)
  qmin = bnp.min(_flatten_cells(g), axis=1)
  qmax = bnp.max(_flatten_cells(g), axis=1)
  c = _flatten_cells(_with_levels(fv_grid["spheremp_fv"], f))
  out = _caas_clip_and_sum(_flatten_cells(f), c, qmin, qmax)
  return bnp.reshape(out, f.shape)


def _with_levels(a, like):
  """Broadcast a per-point array ``(E, a, a)`` against ``like``'s shape."""
  return bnp.expand_dims(a, -1) if like.ndim == 4 else a


# ---------------------------------------------------------------------------
# Vector remap through the element map (gfr_g2f_vector / gfr_f2g_vector)
# ---------------------------------------------------------------------------
def _ref_g2f(field, fv_grid):
  """Unweighted (metric-free) reference-element GLL -> FV average."""
  m = fv_grid["subcell_integral"]                                # (nf, np)
  return bnp.einsum("ci,dj,eij...->ecd...", m, m, field) / fv_grid["w_ff"]


def _ref_f2g(field_fv, fv_grid):
  """Unweighted (metric-free) reference-element FV -> GLL remap."""
  r = fv_grid["fv_to_gll_ref"]                                   # (np, nf)
  return bnp.einsum("ic,jd,ecd...->eij...", r, r, field_fv)


def _apply_map(tensor, u):
  """``out_a = sum_b tensor[..., a, b] u[..., b]`` with the tensor broadcast
  over any level axis of ``u`` (``u``: ``(E, a, a[, K], 2)``)."""
  t = tensor if u.ndim == tensor.ndim - 1 else bnp.expand_dims(tensor, 3)
  return bnp.einsum("...ab,...b->...a", t, u)


def gll_to_fv_vector(u, h_grid, fv_grid):
  """GLL -> FV remap of a physical (lon-lat) vector (``gfr_g2f_vector``).

  The vector is converted to contravariant (reference-element) components
  with ``h_grid["physical_to_contra"]``, each component is remapped with the
  *unweighted* reference operator, and the result is converted back to
  physical components with the FV-point tensors -- so the remap never sees
  the pole-crossing component gradients.  ``u`` is ``(E, np, np[, K], 2)``
  with pySEs' component convention (as in
  ``operations_2d.operators.physical_to_contravariant``).
  """
  w = _apply_map(h_grid["physical_to_contra"], bnp.flip(u, -1))
  w_fv = _ref_g2f(w, fv_grid)
  return bnp.flip(_apply_map(fv_grid["contra_to_physical_fv"], w_fv), -1)


def fv_to_gll_vector(u_fv, h_grid, fv_grid):
  """FV -> GLL remap of a physical vector (``gfr_f2g_vector``); the exact
  mirror of :func:`gll_to_fv_vector`."""
  w_fv = _apply_map(fv_grid["physical_to_contra_fv"], bnp.flip(u_fv, -1))
  w = _ref_f2g(w_fv, fv_grid)
  return bnp.flip(_apply_map(h_grid["contra_to_physical"], w), -1)


# ---------------------------------------------------------------------------
# Hydrostatic FV pressure (calc_dp_fv) and midlevel pressure
# ---------------------------------------------------------------------------
def calc_dp_fv(ps_fv, v_grid):
  """Rebuild FV pressure increments hydrostatically from remapped ``ps``
  (HOMME ``calc_dp_fv``): CAM computes ``pdel`` from the hybrid coefficients
  and ``ps``, so the FV grid must too -- directly remapping ``dp`` disagrees
  numerically.  Same formula as ``mass_coordinate.surface_mass_to_d_mass``.
  """
  return surface_mass_to_d_mass(ps_fv, v_grid)


def _p_mid(dp, ptop):
  """Hydrostatic midlevel pressure from layer increments (HOMME
  ``get_field('p')``): ``p_k = ptop + cumsum(dp)_k - dp_k/2``."""
  return ptop + bnp.cumsum(dp, axis=-1) - 0.5 * dp


# ---------------------------------------------------------------------------
# Full-state drivers (gfr_dyn_to_fv_phys / gfr_fv_phys_to_dyn)
# ---------------------------------------------------------------------------
def dyn_to_fv_phys(ps, phis, T, uv, omega, q, dp_gll, h_grid, fv_grid,
                   v_grid, kappa=KAPPA_DRY):
  """Remap the dynamics state to the FV physics grid (``gfr_dyn_to_fv_phys``).

  Parameters
  ----------
  ps, phis : ``(E, np, np)``; T, omega : ``(E, np, np, K)``;
  uv : ``(E, np, np, K, 2)`` physical winds; q : dict of ``(E, np, np, K)``
  mixing ratios; dp_gll : ``(E, np, np, K)`` layer mass consistent with
  ``ps``; v_grid : vertical grid struct (hybrid coefficients).

  Returns
  -------
  dict with ``"ps"``, ``"phis"`` (limited), ``"dp"`` (hydrostatically
  reconstructed), ``"T"`` (theta-form dp-weighted), ``"uv"`` (map-tensor
  remap), ``"omega"`` (area average) and ``"q"`` (dict, shape-preserving
  dp-weighted remap) on ``(E, nf, nf, ...)``.
  """
  p0 = v_grid["reference_surface_mass"]
  ptop = p0 * v_grid["hybrid_a_i"][0]

  ps_fv = gll_to_fv(ps, fv_grid)
  dp_fv = calc_dp_fv(ps_fv, v_grid)
  phis_fv = gll_to_fv_limited(phis, fv_grid)

  p = _p_mid(dp_gll, ptop)
  p_fv = _p_mid(dp_fv, ptop)
  theta_fv = gll_to_fv_dp(T * (p0 / p) ** kappa, dp_gll, dp_fv, fv_grid)
  T_fv = theta_fv * (p_fv / p0) ** kappa

  return {
      "ps": ps_fv,
      "phis": phis_fv,
      "dp": dp_fv,
      "T": T_fv,
      "uv": gll_to_fv_vector(uv, h_grid, fv_grid),
      "omega": gll_to_fv(omega, fv_grid),
      "q": {name: gll_to_fv_mixing_ratio(qi, dp_gll, dp_fv, fv_grid)
            for name, qi in q.items()},
  }


def fv_phys_to_dyn(dT_fv, duv_fv, q_fv, ps, q_gll, dp_gll, h_grid, fv_grid,
                   v_grid, kappa=KAPPA_DRY, neighbor_minmax=None):
  """Return the physics update to the dynamics grid (``gfr_fv_phys_to_dyn``).

  HOMME's asymmetric convention: temperature and winds come back as
  *tendencies* (increments over the physics step), tracers as *full updated
  states*.

  Parameters
  ----------
  dT_fv : ``(E, nf, nf, K)``
      Physics temperature increment on the FV grid.
  duv_fv : ``(E, nf, nf, K, 2)``
      Physics wind increment (physical components).
  q_fv : dict of ``(E, nf, nf, K)``
      UPDATED tracer states on the FV grid.
  ps, q_gll, dp_gll
      Pre-physics dynamics-grid surface pressure, tracer states (dict) and
      layer mass (the same arrays passed to :func:`dyn_to_fv_phys`).
  neighbor_minmax : callable, optional
      ``(qmin, qmax) -> (qmin, qmax)`` on ``(E, K)`` arrays, expanding the
      per-element FV tracer bounds over element neighbors (HOMME's
      ``gfr_f2g_mixing_ratios_he`` halo exchange; build one from
      ``operations_2d.local_assembly.minmax_scalar``).  Default: bounds stay
      element-local, which is more (not less) restrictive.

  Returns
  -------
  dict with ``"FT"`` ``(E, np, np, K)`` (temperature increment), ``"FM"``
  ``(E, np, np, K, 2)`` (wind increment) and ``"FQ"`` (dict of updated,
  limited GLL tracer states).  All three are element-discontinuous: apply
  DSS afterwards (HOMME's ``gfr_f2g_dss``); FQ's limiter conserves the
  element ``spheremp*dp`` mass, which is what makes the post-DSS field
  conservative.
  """
  p0 = v_grid["reference_surface_mass"]
  ptop = p0 * v_grid["hybrid_a_i"][0]

  ps_fv = gll_to_fv(ps, fv_grid)
  dp_fv = calc_dp_fv(ps_fv, v_grid)

  p = _p_mid(dp_gll, ptop)
  p_fv = _p_mid(dp_fv, ptop)
  dtheta = dT_fv * (p0 / p_fv) ** kappa
  FT = fv_to_gll_dp(dtheta, dp_fv, dp_gll, fv_grid) * (p / p0) ** kappa

  FM = fv_to_gll_vector(duv_fv, h_grid, fv_grid)

  c = _flatten_cells(dp_gll * _with_levels(fv_grid["spheremp"], dp_gll))
  FQ = {}
  for name, q1_fv in q_fv.items():
    q0 = q_gll[name]
    # FV tendency relative to the (limited) forward remap of the GLL state,
    # remapped back and added -- so remap error cancels in the zero-tendency
    # limit.
    q0_fv = gll_to_fv_mixing_ratio(q0, dp_gll, dp_fv, fv_grid)
    dq = fv_to_gll_dp(q1_fv - q0_fv, dp_fv, dp_gll, fv_grid)
    q1 = q0 + dq
    # Bounds: FV state extrema (optionally neighbor-expanded), augmented
    # with the GLL state bounds so a zero tendency returns q0 exactly.
    qmin = bnp.min(_flatten_cells(q1_fv), axis=1)
    qmax = bnp.max(_flatten_cells(q1_fv), axis=1)
    if neighbor_minmax is not None:
      qmin, qmax = neighbor_minmax(qmin, qmax)
    qmin = bnp.minimum(qmin, bnp.min(_flatten_cells(q0), axis=1))
    qmax = bnp.maximum(qmax, bnp.max(_flatten_cells(q0), axis=1))
    out = _caas_clip_and_sum(_flatten_cells(q1), c, qmin, qmax)
    FQ[name] = bnp.reshape(out, q1.shape)

  return {"FT": FT, "FM": FM, "FQ": FQ}
