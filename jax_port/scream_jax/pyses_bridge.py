"""pySEs <-> SCREAM physics bridge.

Couples the validated SCREAMv1 physics suite (scream_jax.driver) to the
pySEs dynamical core, playing the role EAMxx's HommeDynamics interface +
SurfaceCouplingImporter play in SCREAM:

    (t, state) = simulator-yield
    fields   = pyses_to_eamxx(state-derived columns, coupler state)
    fields'  = ScreamPhysics.step(fields, dt, ...)
    forcing  = eamxx_to_pyses_forcing(fields, fields', dt)
    simulator.send(forcing)   # applied via advance_coupling_step

Conventions bridged:
  - pySEs arrays are (elem, np, np, lev) on GLL points with DRY layer
    mass d_mass [Pa] + DRY-based mixing ratios; EAMxx fields are
    (ncol, lev) with WET pseudo_density and WET mixing ratios. pySEs'
    horizontal_wind is already PHYSICAL lon-lat (component 0 = zonal u,
    1 = meridional v, m/s) — its operators convert to contravariant
    internally — so wind components map to EAMxx horiz_winds directly
    with no coordinate transform.
  - wet dp = dry dp * (1 + sum of dry-based water mixing ratios), and
    q_wet = q_dry * dp_dry / dp_wet — the same dp-based conversions the
    EAMxx process pre/post steps use.
  - Physics-internal state EAMxx carries between steps but pySEs does
    not (tke, sgs_buoy_flux, cldfrac_liq, eddy_diff_mom, hydrometeor
    numbers, qv/T_prev, precip masses, rad heating) lives in the
    coupler and persists across calls.

Only mass (hydrostatic, CAM-SE 'T' thermodynamic variable) coupling is
implemented; theta-based HOMME states convert to T first on the pySEs
side.

Physics-grid options:
  - PysesScreamCoupler runs physics directly on the GLL columns (np4
    physics; SCREAM without pg2).
  - PysesScreamCouplerPg2 runs physics on the pg-N FV grid the way
    SCREAM does operationally: the GLL state is remapped with the
    operational HOMME gllfvremap layer
    (pyses_ext.finite_volume_grid_operational — dp-weighted theta-form
    remaps, CAAS-limited tracers, map-tensor winds), physics runs on the
    E*nf*nf FV columns, and the increments return to GLL with HOMME's
    tendency(T,uv)/state(q) asymmetry. Requires pySEs importable.
"""

import numpy as np

from .driver import ScreamPhysics

# water species the wet pressure accounts for (P3's prognostic set)
WATER_SPECIES = ("qv", "qc", "qi", "qr", "qm")
GRAVIT = 9.80616


def columnize(a):
    """(elem, np, np, [lev...]) -> (ncol, [lev...])."""
    a = np.asarray(a)
    return a.reshape((-1,) + a.shape[3:])


def decolumnize(a, elem_shape):
    """(ncol, [lev...]) -> (elem, np, np, [lev...])."""
    a = np.asarray(a)
    return a.reshape(tuple(elem_shape) + a.shape[1:])


def pyses_to_eamxx(T, u, v, omega, d_mass_dry, q_dry, ptop):
    """Build the EAMxx thermodynamic/moisture fields from columnized
    pySEs state.

    T, u, v, omega, d_mass_dry: (ncol, nlev); q_dry: dict of DRY-based
    mixing ratios (must include 'qv'; other WATER_SPECIES default 0);
    ptop: model-top pressure [Pa]. Returns a field dict with p_int,
    p_mid, p_dry_int, p_dry_mid, pseudo_density(_dry), T_mid,
    horiz_winds, omega and WET mixing ratios for every q_dry entry."""
    d_mass_dry = np.asarray(d_mass_dry, dtype=np.float64)
    ncol, nlev = d_mass_dry.shape

    water = np.zeros((ncol, nlev))
    for s in WATER_SPECIES:
        if s in q_dry:
            water = water + np.asarray(q_dry[s], dtype=np.float64)
    dp_wet = d_mass_dry * (1.0 + water)

    def to_int(dp):
        p_int = np.zeros((ncol, nlev + 1))
        p_int[:, 0] = ptop
        p_int[:, 1:] = ptop + np.cumsum(dp, axis=1)
        return p_int

    p_dry_int = to_int(d_mass_dry)
    p_int = to_int(dp_wet)

    out = {
        "pseudo_density": dp_wet,
        "pseudo_density_dry": d_mass_dry,
        "p_int": p_int,
        "p_mid": 0.5 * (p_int[:, :-1] + p_int[:, 1:]),
        "p_dry_int": p_dry_int,
        "p_dry_mid": 0.5 * (p_dry_int[:, :-1] + p_dry_int[:, 1:]),
        "T_mid": np.asarray(T, dtype=np.float64).copy(),
        "omega": np.asarray(omega, dtype=np.float64).copy(),
    }
    hw = np.zeros((ncol, 2, nlev))
    hw[:, 0, :] = np.asarray(u, dtype=np.float64)
    hw[:, 1, :] = np.asarray(v, dtype=np.float64)
    out["horiz_winds"] = hw

    for name, q in q_dry.items():
        out[name] = np.asarray(q, dtype=np.float64) * d_mass_dry / dp_wet
    return out


def eamxx_to_pyses_forcing(before, after, q_dry_before, dt):
    """Physics tendencies in pySEs variables.

    before/after: EAMxx field dicts around the physics step;
    q_dry_before: the pre-physics DRY mixing ratios (dict).
    Returns dict with 'FT' (K/s), 'FU'/'FV' (m/s^2) and 'FQ' (dict of
    DRY mixing-ratio tendencies, 1/s). Physics conserves dry mass, so
    dry->wet conversion uses the (fixed) dry dp with the UPDATED wet dp
    implied by the new water content."""
    dp_dry = np.asarray(before["pseudo_density_dry"])

    # new wet dp implied by post-physics wet mixing ratios:
    # q_wet = q_dry*dp_dry/dp_wet and dp_wet = dp_dry + sum(q_dry*dp_dry)
    # => dp_wet_new = dp_dry / (1 - sum(q_wet_new))
    total_q_wet = np.zeros_like(dp_dry)
    for s in WATER_SPECIES:
        if s in after:
            total_q_wet = total_q_wet + np.asarray(after[s])
    dp_wet_new = dp_dry / (1.0 - total_q_wet)

    FQ = {}
    for name, q0 in q_dry_before.items():
        q_dry_new = np.asarray(after[name]) * dp_wet_new / dp_dry
        FQ[name] = (q_dry_new - np.asarray(q0)) / dt

    FT = (np.asarray(after["T_mid"]) - np.asarray(before["T_mid"])) / dt
    FU = (np.asarray(after["horiz_winds"])[:, 0, :]
          - np.asarray(before["horiz_winds"])[:, 0, :]) / dt
    FV = (np.asarray(after["horiz_winds"])[:, 1, :]
          - np.asarray(before["horiz_winds"])[:, 1, :]) / dt
    return {"FT": FT, "FU": FU, "FV": FV, "FQ": FQ}


class PysesScreamCoupler:
    """Persistent physics state + the suite driver, for use in a pySEs
    time loop."""

    #: physics-internal prognostic/persistent fields with their initial
    #: values (EAMxx initializes these from the IC file or to constants)
    INTERNAL_2D = {
        "tke": 0.0004, "sgs_buoy_flux": 0.0, "eddy_diff_mom": 0.0,
        "cldfrac_liq": 0.0, "nc": 0.0, "nr": 0.0, "ni": 0.0, "bm": 0.0,
        "qm": 0.0, "rad_heating_pdel": 0.0, "nc_nuceat_tend": 0.0,
        "ni_activated": 0.0, "inv_qc_relvar": 1.0,
    }
    INTERNAL_1D = {"precip_liq_surf_mass": 0.0, "precip_ice_surf_mass": 0.0}

    def __init__(self, data_dir, hyam, hybm, lat_deg, lon_deg, cell_length,
                 ncol, nlev, surface, o3_vmr, phis=None,
                 mac_mic_subcycles=6, year=2021, **phys_kw):
        """surface: dict with surf_evap, surf_sens_flux, surf_mom_flux
        (ncol,2), surf_lw_flux_up, sfc_alb_{dir,dif}_{vis,nir} (ncol,) —
        prescribed or supplied by a surface model each step.
        o3_vmr: (ncol, nlev) prescribed ozone."""
        self.phys = ScreamPhysics(data_dir, hyam, hybm, lat_deg, lon_deg,
                                  cell_length,
                                  mac_mic_subcycles=mac_mic_subcycles,
                                  year=year, **phys_kw)
        self.state = {}
        for k, v in self.INTERNAL_2D.items():
            self.state[k] = np.full((ncol, nlev), v)
        for k, v in self.INTERNAL_1D.items():
            self.state[k] = np.full(ncol, v)
        self.state["o3_volume_mix_ratio"] = np.asarray(o3_vmr,
                                                       dtype=np.float64)
        self.state["phis"] = (np.zeros(ncol) if phis is None
                              else np.asarray(phis))
        self.surface = {k: np.asarray(v, dtype=np.float64)
                        for k, v in surface.items()}
        self._first = True

    def step(self, T, u, v, omega, d_mass_dry, q_dry, ptop, dt, nstep,
             doy_start, surface=None):
        """One physics step from columnized pySEs prognostics. Returns
        (forcing dict for advance_coupling_step's physics_forcing,
        diagnostics dict of all EAMxx fields)."""
        if surface is not None:
            self.surface = {k: np.asarray(v, dtype=np.float64)
                            for k, v in surface.items()}

        # coupler-internal state first, then pySEs-derived fields (so any
        # species pySEs advects — including hydrometeor numbers if present
        # in q_dry — override the internally persisted copies), then the
        # surface data
        fields = dict(self.state)
        fields.update(pyses_to_eamxx(T, u, v, omega, d_mass_dry, q_dry,
                                     ptop))
        fields.update(self.surface)
        if self._first:
            # qv/T at the previous micro step start equal to current
            fields.setdefault("qv_prev_micro_step", fields["qv"].copy())
            fields.setdefault("T_prev_micro_step", fields["T_mid"].copy())
            self._first = False

        before = {k: np.array(v) for k, v in fields.items()
                  if k in ("T_mid", "horiz_winds", "pseudo_density_dry")}
        out = self.phys.step(fields, dt, nstep, doy_start)

        # persist physics-internal state for the next call
        for k in list(self.INTERNAL_2D) + list(self.INTERNAL_1D) + [
                "qv_prev_micro_step", "T_prev_micro_step"]:
            if k in out:
                self.state[k] = np.asarray(out[k])

        forcing = eamxx_to_pyses_forcing(before, out, q_dry, dt)
        return forcing, out


class PysesScreamCouplerPg2:
    """pg-N (SCREAM-operational) variant of :class:`PysesScreamCoupler`.

    Physics runs on the element-local FV grid (``E*nf*nf`` columns)
    instead of the GLL points, with the state remapped both ways by the
    operational gllfvremap layer.  The remap works on pySEs' DRY state
    (dry layer mass follows the hybrid coordinate, which is what the
    hydrostatic dp reconstruction assumes); the dry->wet conversion the
    EAMxx physics needs happens per-column on the FV grid, inside the
    inner coupler.  SCREAM remaps the wet state instead — an accepted
    divergence that keeps the bridge consistent with pySEs' dry-mass
    prognostics.

    Returned forcings are GLL-shaped (E, np, np, lev) tendencies for
    advance_coupling_step, element-discontinuous like HOMME's before its
    DSS pass: supply ``dss`` (a callable applied to each forcing field)
    or project on the pySEs side.
    """

    def __init__(self, data_dir, h_grid, fv_grid, v_grid, lat_deg, lon_deg,
                 cell_length, surface, o3_vmr, phis_gll=None, phis_fv=None,
                 dss=None, neighbor_minmax=None, **coupler_kw):
        """h_grid/fv_grid/v_grid: pySEs horizontal grid (with the
        contra/physical map tensors), extended FV grid
        (init_fv_grid + extend_fv_grid_operational) and vertical grid.
        lat_deg/lon_deg/cell_length/surface/o3_vmr: per FV COLUMN
        (E*nf*nf, row-major over (elem, nf, nf)).
        phis_fv: (E, nf, nf) FV topography (SCREAM's preferred source);
        phis_gll: (E, np, np) alternative, remapped with the limited
        scalar remap. dss: optional callable field -> field applied to
        every returned forcing. neighbor_minmax: optional (qmin, qmax)
        -> (qmin, qmax) halo-exchange hook for the tracer limiter."""
        from pyses_ext import finite_volume_grid_operational as opfv
        self._opfv = opfv
        self.h_grid, self.fv_grid, self.v_grid = h_grid, fv_grid, v_grid
        self.nf = int(fv_grid["nf"])
        self.npt = int(fv_grid["npt"])
        self.num_elem = int(fv_grid["num_elem"])
        self.elem_shape = (self.num_elem, self.nf, self.nf)
        ncol = self.num_elem * self.nf * self.nf

        hyam = np.asarray(v_grid["hybrid_a_m"], dtype=np.float64)
        hybm = np.asarray(v_grid["hybrid_b_m"], dtype=np.float64)
        nlev = hyam.size
        p0 = float(v_grid["reference_surface_mass"])
        self.ptop = p0 * float(np.asarray(v_grid["hybrid_a_i"])[0])
        self._zero_phis_gll = np.zeros((self.num_elem, self.npt, self.npt))

        if phis_fv is not None:
            phis_cols = np.asarray(phis_fv, dtype=np.float64).reshape(ncol)
        elif phis_gll is not None:
            phis_cols = columnize(np.asarray(opfv.gll_to_fv_limited(
                np.asarray(phis_gll, dtype=np.float64), fv_grid)))
        else:
            phis_cols = np.zeros(ncol)

        self.inner = PysesScreamCoupler(
            data_dir, hyam, hybm, lat_deg, lon_deg, cell_length,
            ncol, nlev, surface, o3_vmr, phis=phis_cols, **coupler_kw)
        self.dss = dss
        self.neighbor_minmax = neighbor_minmax

    def step(self, T, horizontal_wind, omega, d_mass_dry, q_dry, dt, nstep,
             doy_start, surface=None):
        """One physics step from ELEMENT-SHAPED pySEs prognostics:
        T/omega/d_mass_dry (E, np, np, lev), horizontal_wind
        (E, np, np, lev, 2) physical (u, v), q_dry dict of DRY mixing
        ratios (E, np, np, lev). Returns (forcing dict with GLL-shaped
        FT [K/s], FU/FV [m/s^2], FQ [dict, 1/s], all pre-DSS unless
        ``dss`` was supplied; FV-column diagnostics dict)."""
        opfv = self._opfv
        T = np.asarray(T, dtype=np.float64)
        hw = np.asarray(horizontal_wind, dtype=np.float64)
        omega = np.asarray(omega, dtype=np.float64)
        dp_dry = np.asarray(d_mass_dry, dtype=np.float64)
        q_dry = {k: np.asarray(v, dtype=np.float64)
                 for k, v in q_dry.items()}
        ps_dry = self.ptop + dp_dry.sum(axis=-1)

        # ---- dynamics -> FV physics columns (dry state) ----
        fwd = opfv.dyn_to_fv_phys(ps_dry, self._zero_phis_gll, T, hw, omega,
                                  q_dry, dp_dry, self.h_grid, self.fv_grid,
                                  self.v_grid)
        uv_fv = np.asarray(fwd["uv"])
        q_fv_cols = {k: columnize(np.asarray(v))
                     for k, v in fwd["q"].items()}

        forcing_cols, out = self.inner.step(
            columnize(np.asarray(fwd["T"])), columnize(uv_fv[..., 0]),
            columnize(uv_fv[..., 1]), columnize(np.asarray(fwd["omega"])),
            columnize(np.asarray(fwd["dp"])), q_fv_cols, self.ptop,
            dt, nstep, doy_start, surface)

        # ---- FV physics increments -> dynamics grid ----
        dT_fv = decolumnize(dt * forcing_cols["FT"], self.elem_shape)
        duv_fv = decolumnize(
            np.stack([dt * forcing_cols["FU"], dt * forcing_cols["FV"]],
                     axis=-1), self.elem_shape)
        q1_fv = {k: decolumnize(q_fv_cols[k] + dt * forcing_cols["FQ"][k],
                                self.elem_shape)
                 for k in q_fv_cols}
        back = opfv.fv_phys_to_dyn(dT_fv, duv_fv, q1_fv, ps_dry, q_dry,
                                   dp_dry, self.h_grid, self.fv_grid,
                                   self.v_grid,
                                   neighbor_minmax=self.neighbor_minmax)

        FM = np.asarray(back["FM"])
        forcing = {
            "FT": np.asarray(back["FT"]) / dt,
            "FU": FM[..., 0] / dt,
            "FV": FM[..., 1] / dt,
            "FQ": {k: (np.asarray(back["FQ"][k]) - q_dry[k]) / dt
                   for k in q1_fv},
        }
        if self.dss is not None:
            forcing = {
                "FT": self.dss(forcing["FT"]),
                "FU": self.dss(forcing["FU"]),
                "FV": self.dss(forcing["FV"]),
                "FQ": {k: self.dss(v) for k, v in forcing["FQ"].items()},
            }
        return forcing, out
