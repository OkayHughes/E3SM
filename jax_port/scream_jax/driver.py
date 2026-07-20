"""SCREAMv1 physics-suite driver: the validated process chain assembled
in EAMxx AtmosphereDriver order.

Mirrors the physics-only coupled configuration
(components/eamxx/tests/multi-process/physics_only/shoc_cld_spa_p3_rrtmgp):

    atm_procs_list: [mac_mic, rrtmgp]
    mac_mic (subcycled N times at dt/N): shoc -> cld_fraction -> spa -> p3
    rrtmgp (once per step)

Each sub-step is the exact swap-tested process step from
scream_jax/{shoc,cld_fraction,spa,p3,rrtmgp}/. State flows through a
plain dict of EAMxx field names (the same fields the EAMxx field
manager holds), so this driver is directly comparable against a
multi-process pyeamxx golden run — and is the unit the pySEs bridge
(scream_jax/pyses_bridge.py) drives.

Time bookkeeping matches EAMxx: subcycled processes see dt/N and
timestamps advancing by dt/N; SPA uses the END-of-substep timestamp;
RRTMGP uses the START-of-step timestamp and step counter.
"""

import numpy as np

from .shoc import shoc_init
from .shoc.process import shoc_process_step
from .cld_fraction.main import cld_fraction_main
from .spa.process import spa_process_step, load_spa_data
from .p3 import DEFAULT_OPTS as P3_DEFAULT_OPTS
from .p3 import tables as p3_tables
from .p3.process import p3_process_step
from .rrtmgp.coefficients import load_kdist
from .rrtmgp.cloud_optics import load_cloud_optics
from .rrtmgp.process import (DEFAULT_PARAMS as RRTMGP_DEFAULT_PARAMS,
                             GAS_NAMES, rrtmgp_process_step)

P0 = 100000.0

SHOC_PARAMS = dict(lambda_low=0.001, lambda_high=0.08, lambda_slope=2.65,
                   lambda_thresh=0.02, thl2tune=1.0, qw2tune=1.0,
                   qwthl2tune=1.0, w2tune=1.0, length_fac=0.5,
                   c_diag_3rd_mom=7.0, ckh=0.1, ckm=0.1)
CLDFRAC_ICE_THRESHOLD = 1e-12
CLDFRAC_ICE_4OUT_THRESHOLD = 1e-5


class ScreamPhysics:
    """Holds static data (tables, grids, params) and steps the suite."""

    def __init__(self, data_dir, hyam, hybm, lat_deg, lon_deg, cell_length,
                 mac_mic_subcycles=6, p3_opts=None, rrtmgp_params=None,
                 year=2021, spa_col_indices=None):
        """data_dir: directory with the p3 tables/, rrtmgp coefficient
        files and the SPA data file (e3sm-inputdata/atm/scream layout).
        hyam/hybm: hybrid midpoint coefficients (for SHOC's npbl).
        cell_length: (ncol,) grid spacing for SHOC [m].
        spa_col_indices: optional (ncol,) indices selecting, for each
        physics column, its source column in the SPA data file — for
        running on a grid other than the file's (e.g. pg2 physics with
        an np4 SPA file; operationally SCREAM ships a per-grid file)."""
        d = str(data_dir)
        self.p3_tables = p3_tables.p3_init(f"{d}/tables")
        self.p3_opts = dict(P3_DEFAULT_OPTS, **(p3_opts or {}))
        self.kd_sw = load_kdist(f"{d}/init/rrtmgp-data-sw-g112-210809.nc",
                                GAS_NAMES)
        self.kd_lw = load_kdist(f"{d}/init/rrtmgp-data-lw-g128-210809.nc",
                                GAS_NAMES)
        self.co_sw = load_cloud_optics(
            f"{d}/init/rrtmgp-cloud-optics-coeffs-sw.nc")
        self.co_lw = load_cloud_optics(
            f"{d}/init/rrtmgp-cloud-optics-coeffs-lw.nc")
        self.spa_data = load_spa_data(
            f"{d}/init/spa_file_unified_and_complete_ne2np4L72_20231222.nc")
        if spa_col_indices is not None:
            idx = np.asarray(spa_col_indices)
            ncol_data = self.spa_data["PS"].shape[1]
            self.spa_data = {
                k: (v[:, idx] if (getattr(v, "ndim", 0) >= 2
                                  and v.shape[1] == ncol_data) else v)
                for k, v in self.spa_data.items()}
        self.rrtmgp_params = dict(RRTMGP_DEFAULT_PARAMS,
                                  orbital_year=1990,
                                  **(rrtmgp_params or {}))
        nlev = len(hyam)
        self.npbl = shoc_init(nlev, 0, P0 * (np.asarray(hyam)
                                             + np.asarray(hybm)))
        self.lat = np.asarray(lat_deg)
        self.lon = np.asarray(lon_deg)
        self.cell_length = np.asarray(cell_length)
        self.nsub = mac_mic_subcycles
        self.year = year

    def step(self, s, dt, nstep, doy_start, ad_mode=False,
             spa_outs=None, ad_checkpoint="process"):
        """Advance the physics suite one step.

        s: dict of EAMxx fields (updated in place semantically — a new
        dict is returned). doy_start: 0-based fractional day-of-year of
        the START of the step. nstep: completed-step count at start.
        Returns the updated field dict.

        ad_mode=False (default) is the validated, bitwise-sacred path
        (host numpy round-trips between processes). ad_mode=True runs
        the same physics as one end-to-end traceable JAX computation
        (see _step_ad) so the step can sit under jax.grad / jax.jvp /
        jax.jit; spa_outs / ad_checkpoint are only used in that mode."""
        if ad_mode:
            return self._step_ad(s, dt, nstep, doy_start,
                                 spa_outs=spa_outs,
                                 checkpoint=ad_checkpoint)
        s = dict(s)
        dt_sub = dt / self.nsub

        for isub in range(self.nsub):
            # ---- shoc ----
            out = shoc_process_step(
                dt_sub, self.npbl, self.cell_length,
                SHOC_PARAMS["lambda_low"], SHOC_PARAMS["lambda_high"],
                SHOC_PARAMS["lambda_slope"], SHOC_PARAMS["lambda_thresh"],
                SHOC_PARAMS["thl2tune"], SHOC_PARAMS["qw2tune"],
                SHOC_PARAMS["qwthl2tune"], SHOC_PARAMS["w2tune"],
                SHOC_PARAMS["length_fac"], SHOC_PARAMS["c_diag_3rd_mom"],
                SHOC_PARAMS["ckh"], SHOC_PARAMS["ckm"], False, False,
                s["T_mid"], s["p_mid"], s["p_int"], s["pseudo_density"],
                s["omega"], s["phis"], s["surf_sens_flux"], s["surf_evap"],
                s["surf_mom_flux"][:, 0], s["surf_mom_flux"][:, 1],
                s["qv"], s["qc"], s["tke"],
                s["horiz_winds"][:, 0, :], s["horiz_winds"][:, 1, :],
                s["cldfrac_liq"], s["sgs_buoy_flux"], s["eddy_diff_mom"])
            for k, v in out.items():
                if k in ("u_wind", "v_wind"):
                    continue
                s[k] = np.asarray(v)
            hw = np.array(s["horiz_winds"])
            hw[:, 0, :] = np.asarray(out["u_wind"])
            hw[:, 1, :] = np.asarray(out["v_wind"])
            s["horiz_winds"] = hw

            # ---- cld_fraction ----
            ice, tot, ice4, tot4 = cld_fraction_main(
                CLDFRAC_ICE_THRESHOLD, CLDFRAC_ICE_4OUT_THRESHOLD,
                s["qi"], s["cldfrac_liq"])
            s["cldfrac_ice"] = np.asarray(ice)
            s["cldfrac_tot"] = np.asarray(tot)
            s["cldfrac_ice_for_analysis"] = np.asarray(ice4)
            s["cldfrac_tot_for_analysis"] = np.asarray(tot4)

            # ---- spa (end-of-substep timestamp) ----
            doy_end_sub = doy_start + (isub + 1) * dt_sub / 86400.0
            spa_out = spa_process_step(self.spa_data, doy_end_sub, s["p_mid"])
            for k, v in spa_out.items():
                s[k] = np.asarray(v)

            # ---- p3 (prescribed CCN from spa's nccn) ----
            p3_out = p3_process_step(
                dt_sub, True, True, True, False, False, False, False, False,
                s["T_mid"], s["p_mid"], s["p_dry_mid"],
                s["pseudo_density"], s["pseudo_density_dry"],
                s["cldfrac_tot"],
                s["qv"], s["qc"], s["nc"], s["qr"], s["nr"],
                s["qi"], s["qm"], s["ni"], s["bm"],
                s["qv_prev_micro_step"], s["T_prev_micro_step"],
                s["nc_nuceat_tend"], s["nccn"], s["ni_activated"],
                s["inv_qc_relvar"],
                s["precip_liq_surf_mass"], s["precip_ice_surf_mass"],
                self.p3_tables, self.p3_opts)
            for k, v in p3_out.items():
                s[k] = np.asarray(v)

        # ---- rrtmgp (start-of-step timestamp) ----
        rad_out = rrtmgp_process_step(
            self.kd_sw, self.kd_lw, self.co_sw, self.co_lw,
            self.rrtmgp_params, dt, nstep, self.year, doy_start + 1,
            self.lat, self.lon,
            s["T_mid"], s["p_mid"], s["p_int"], s["pseudo_density"],
            s["sfc_alb_dir_vis"], s["sfc_alb_dir_nir"],
            s["sfc_alb_dif_vis"], s["sfc_alb_dif_nir"],
            s["qv"], s["qc"], s["nc"], s["qi"], s["cldfrac_tot"],
            s["eff_radius_qc"], s["eff_radius_qi"], s["surf_lw_flux_up"],
            s["o3_volume_mix_ratio"], s["rad_heating_pdel"],
            # do_aerosol_rad (default true): aerosol optics from SPA
            aero_tau_sw=s["aero_tau_sw"], aero_ssa_sw=s["aero_ssa_sw"],
            aero_g_sw=s["aero_g_sw"], aero_tau_lw=s["aero_tau_lw"])
        for k, v in rad_out.items():
            s[k] = np.asarray(v)
        return s

    # ------------------------------------------------------------------
    # Differentiable (traceable) path
    # ------------------------------------------------------------------
    def spa_substep_outputs(self, dt, doy_start, p_mid):
        """SPA outputs for every mac_mic substep of one step, computed
        host-side (numpy). SPA is pure data interpolation in TIME and
        source->target PRESSURE; no process in the suite modifies p_mid,
        so hoisting it out of the traced region is exact. p_mid must be
        concrete (numpy or non-traced jax array)."""
        dt_sub = dt / self.nsub
        p = np.asarray(p_mid)
        return [spa_process_step(self.spa_data,
                                 doy_start + (isub + 1) * dt_sub / 86400.0,
                                 p)
                for isub in range(self.nsub)]

    def _step_ad(self, s, dt, nstep, doy_start, spa_outs=None,
                 checkpoint="process"):
        """One suite step as a single traceable JAX computation.

        Same physics as the default path, with exactly these changes:
          (a) no inter-process host round-trips — state stays jax
              arrays end-to-end, so jax.grad/jvp/jit work through the
              whole step;
          (b) SPA is hoisted host-side (spa_substep_outputs; exact,
              see there). Pass spa_outs precomputed if s["p_mid"] is a
              tracer at call time;
          (c) P3 sedimentation uses the bounded-scan realization
              (sed_use_while_loop=False; bitwise-equal primal,
              reverse-mode differentiable);
          (d) process calls are wrapped in jax.checkpoint so reverse
              mode stores only process-boundary states.

        checkpoint: "process" (each of shoc/cld_fraction/p3/rrtmgp; the
        default), "subcycle" (one checkpoint per mac_mic substep body +
        one for rrtmgp; fewer residuals, more recompute), or "none".
        RRTMGP's orbital/zenith pieces stay host-side inside the traced
        region: they depend only on time and concrete lat/lon, never on
        state, so they trace as constants."""
        import jax
        import jax.numpy as jnp

        s = dict(s)
        dt_sub = dt / self.nsub
        if spa_outs is None:
            spa_outs = self.spa_substep_outputs(dt, doy_start, s["p_mid"])

        ck = jax.checkpoint if checkpoint != "none" else (lambda f: f)
        ck_proc = jax.checkpoint if checkpoint == "process" else (lambda f: f)
        ck_sub = jax.checkpoint if checkpoint == "subcycle" else (lambda f: f)

        def shoc_fn(sub):
            out = shoc_process_step(
                dt_sub, self.npbl, self.cell_length,
                SHOC_PARAMS["lambda_low"], SHOC_PARAMS["lambda_high"],
                SHOC_PARAMS["lambda_slope"], SHOC_PARAMS["lambda_thresh"],
                SHOC_PARAMS["thl2tune"], SHOC_PARAMS["qw2tune"],
                SHOC_PARAMS["qwthl2tune"], SHOC_PARAMS["w2tune"],
                SHOC_PARAMS["length_fac"], SHOC_PARAMS["c_diag_3rd_mom"],
                SHOC_PARAMS["ckh"], SHOC_PARAMS["ckm"], False, False,
                sub["T_mid"], sub["p_mid"], sub["p_int"],
                sub["pseudo_density"], sub["omega"], sub["phis"],
                sub["surf_sens_flux"], sub["surf_evap"],
                sub["surf_mom_flux"][:, 0], sub["surf_mom_flux"][:, 1],
                sub["qv"], sub["qc"], sub["tke"],
                sub["horiz_winds"][:, 0, :], sub["horiz_winds"][:, 1, :],
                sub["cldfrac_liq"], sub["sgs_buoy_flux"],
                sub["eddy_diff_mom"])
            new = dict(sub)
            for k, v in out.items():
                if k not in ("u_wind", "v_wind"):
                    new[k] = v
            new["horiz_winds"] = jnp.stack(
                [jnp.asarray(out["u_wind"]), jnp.asarray(out["v_wind"])],
                axis=1)
            return new

        def cldfrac_fn(sub):
            ice, tot, ice4, tot4 = cld_fraction_main(
                CLDFRAC_ICE_THRESHOLD, CLDFRAC_ICE_4OUT_THRESHOLD,
                sub["qi"], sub["cldfrac_liq"])
            return dict(sub, cldfrac_ice=ice, cldfrac_tot=tot,
                        cldfrac_ice_for_analysis=ice4,
                        cldfrac_tot_for_analysis=tot4)

        def p3_fn(sub):
            out = p3_process_step(
                dt_sub, True, True, True, False, False, False, False,
                False,
                sub["T_mid"], sub["p_mid"], sub["p_dry_mid"],
                sub["pseudo_density"], sub["pseudo_density_dry"],
                sub["cldfrac_tot"],
                sub["qv"], sub["qc"], sub["nc"], sub["qr"], sub["nr"],
                sub["qi"], sub["qm"], sub["ni"], sub["bm"],
                sub["qv_prev_micro_step"], sub["T_prev_micro_step"],
                sub["nc_nuceat_tend"], sub["nccn"], sub["ni_activated"],
                sub["inv_qc_relvar"],
                sub["precip_liq_surf_mass"], sub["precip_ice_surf_mass"],
                self.p3_tables, self.p3_opts,
                sed_use_while_loop=False)
            new = dict(sub)
            new.update(out)
            return new

        def rad_fn(sub):
            out = rrtmgp_process_step(
                self.kd_sw, self.kd_lw, self.co_sw, self.co_lw,
                self.rrtmgp_params, dt, nstep, self.year, doy_start + 1,
                self.lat, self.lon,
                sub["T_mid"], sub["p_mid"], sub["p_int"],
                sub["pseudo_density"],
                sub["sfc_alb_dir_vis"], sub["sfc_alb_dir_nir"],
                sub["sfc_alb_dif_vis"], sub["sfc_alb_dif_nir"],
                sub["qv"], sub["qc"], sub["nc"], sub["qi"],
                sub["cldfrac_tot"],
                sub["eff_radius_qc"], sub["eff_radius_qi"],
                sub["surf_lw_flux_up"],
                sub["o3_volume_mix_ratio"], sub["rad_heating_pdel"],
                aero_tau_sw=sub["aero_tau_sw"],
                aero_ssa_sw=sub["aero_ssa_sw"],
                aero_g_sw=sub["aero_g_sw"],
                aero_tau_lw=sub["aero_tau_lw"])
            new = dict(sub)
            new.update({k: jnp.asarray(v) for k, v in out.items()})
            return new

        def substep(sub, isub):
            sub = ck_proc(shoc_fn)(sub)
            sub = ck_proc(cldfrac_fn)(sub)
            sub = dict(sub)
            sub.update(spa_outs[isub])
            return ck_proc(p3_fn)(sub)

        for isub in range(self.nsub):
            s = ck_sub(lambda sub, i=isub: substep(sub, i))(s)
        return ck(rad_fn)(s)
