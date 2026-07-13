#!/usr/bin/env python3
"""Generate CLUBB goldens (run in the scream-dev container).

Drives the REAL CLUBB Fortran (eam_clubb_f, the full unmodified stack,
see build_clubb.py) exactly as EAM does and records golden archives:

  golden/clubb_grid.npz         setup_grid arrays + zt2zm/zm2zt/ddzt/
                                ddzm on regime-sweeping test fields,
                                on 3 grids (EAM-like 73-level, coarse
                                41-level, irregular).
  golden/clubb_tridag.npz       tridag_solve (LAPACK dgtsv) systems:
                                diagonally-dominant CLUBB-like,
                                pivot-exercising non-dominant, and
                                multi-RHS cases.
  golden/clubb_sat.npz          sat_mixrat_liq/ice over a (p,T) grid
                                (150-330 K incl. both polynomial-clip
                                branches and the esat~p guard).
  golden/clubb_pdf_closure.npz  pdf_closure (ADG1) over 40 regime-
                                sweeping columns x 73 levels: stable /
                                convective-cloudy / extreme-skewness /
                                degenerate-wp2 / degenerate-xp2 /
                                saturated / cold-cirrus / randomized;
                                plus the Skx and sigma_sqd_w helper
                                inputs+outputs used to build them.
  golden/clubb_pdf_driver.npz   pdf_closure_driver (slice B: the zt+zm
                                double pdf_closure call plumbing,
                                trapezoidal-rule vertical averaging,
                                compute_cloud_cover, clip_rcm) over 40
                                regime-sweeping columns with inputs on
                                their native grids (moments on zm,
                                means on zt), incl. a sharp-dry-notch
                                family that exercises clip_rcm and the
                                cloud-top/base branches of
                                compute_cloud_cover; PLUS 8 "adv"
                                cases whose inputs are REAL
                                advance_clubb_core-advanced states.
                                The verbatim extraction is validated
                                during generation: the 8 adv cases
                                assert BITWISE identity between the
                                extracted pdf_closure_driver replayed
                                on the post-advance state and the
                                passthrough outputs of the real
                                advance_clubb_core (EAMv3
                                ipdf_call_placement=2 places the pdf
                                call last, so its outputs reach the
                                caller untouched).

  golden/clubb_xp2.npz          (slice C) the REAL public
                                advance_xp2_xpyp driven directly (40
                                cases: 16 whose prognostic moments are
                                REAL advance_clubb_core-advanced states
                                + 24 synthetic stress columns hitting
                                every clip: rt_tol^2/thl_tol^2 floors,
                                the 0.5*rtm^2 rtp2 cap, the 1000 up2/
                                vp2 cap, Cauchy-Schwarz rtpthlp
                                clipping, hole-filling via strongly
                                negative forcings, both upwind-ta signs
                                via wp3_on_wp2 sign structure), plus
                                clip_covars_denom cases (post-xp2
                                covariance clipping, l_tke_aniso=T)
                                and the nu2/nu9_vert_res_dep profiles
                                + model-flag configuration.

  golden/clubb_lscale.npz       (slice D) the Lscale/tau
                                infrastructure: the REAL public
                                compute_mixing_length (parcel Lscale,
                                40 columns: physically consistent +
                                stable/unstable/saturated/bitwise-
                                uniform-dCAPE/near-zero-TKE/mu-sweep
                                stress families), the REAL
                                calc_brunt_vaisala_freq_sqd (all three
                                formula variants via flag toggling),
                                calc_stability_correction,
                                term_wp2_splat/term_wp3_splat (C=0
                                EAMv3 + nonzero clip sweep),
                                calc_surface_varnce (surface-flux
                                sign/magnitude sweeps incl. both
                                splat-correction arms), and the full
                                drv_lscale_tau_segment (a verbatim
                                transcription of advance_clubb_core's
                                inline em->thvm->Lscale->tau->Kh->
                                splat->surface->tau_N2 segment,
                                validated BITWISE against the real
                                advance_clubb_core khzm/khzt
                                diagnostics at generation).

Run `python3 gen_clubb_golden.py [archive ...]` with archive names
(grid/tridag/sat/pdf_closure/pdf_driver/xp2/lscale) to regenerate a
subset; no arguments regenerates everything.

EAMv3 tunable parameters: CLUBB compiled-in defaults from the real
read_parameters(-99) with the namelist_defaults_eam.xml phys="default"
clubb_param_nl overrides applied, INCLUDING the coupling rules
read_parameters applies when a value is namelist-set (verbatim from
parameters_tunable.F90): clubb_C2rt set -> C2thl = C2rt and
C2rtthl = C2rt * 1.3 (clubb_c2thl/clubb_c2rtthl are NOT in the EAMv3
namelist); clubb_C6rt NOT set -> C6rt/C6thl keep their defaults.
The full packed params vector is recorded in the golden metadata.
"""
import json
import subprocess
from pathlib import Path

import numpy as np

from fbuild import REPO, import_extension

GOLDEN = Path(__file__).resolve().parent.parent / "golden"
GOLDEN.mkdir(exist_ok=True)

f = import_extension("eam_clubb_f")
d = f.clubb_driver

IDX_NAMES = ["C1", "C1b", "C1c", "C2rt", "C2thl", "C2rtthl", "C6rt",
             "C6rtb", "C6rtc", "C6thl", "C6thlb", "C6thlc", "C7", "C7b",
             "C8", "C11", "C11b", "C11c", "C14", "beta", "gamma_coef",
             "gamma_coefb", "gamma_coefc", "mu", "nu1", "c_K10",
             "c_K10h", "wpxp_L_thresh", "altitude_threshold",
             "Skw_denom_coef", "Skw_max_mag"]

# EAMv3 phys="default" clubb_param_nl values
# (bld/namelist_files/namelist_defaults_eam.xml).  clubb_C2rt=1.75
# triggers the read_parameters coupling C2thl=C2rt, C2rtthl=1.3*C2rt.
# CORRECTIONS in slice C (none affects the outputs of the earlier
# archives -- pdf_closure/pdf_closure_driver read none of these):
#   - C14 = 2.5 (namelist_defaults line ~1968, phys="default"; the
#     CLUBB compiled-in default is 1.0).  advance_xp2_xpyp uses C14
#     in the up2/vp2 dp1/pr1 terms, so slice C REQUIRES it.
#   - c_K10 = c_K10h = 0.35 (lines ~2000-2001; defaults 0.6/1.0);
#     first used by the Kh computation (slice D).
#   - wpxp_L_thresh = 100.0, NOT the generic 60.0: build-namelist's
#     defaults lookup is case-INsensitive (Build::NamelistDefaults
#     get_value lc()'s the name), so the <clubb_wpxp_l_thresh
#     phys="default"> 100.0 entry at line ~2006 beats the
#     attribute-less 60.0 at line 871.  First used by advance_xm_wpxp
#     (slice E).
EAMV3_OVERRIDES = {
    "C1": 2.4, "C1b": 2.8, "C1c": 0.75,
    "C2rt": 1.75, "C2thl": 1.75, "C2rtthl": 1.75 * 1.3,
    "C6rtb": 7.50, "C6rtc": 0.50, "C6thlb": 7.50, "C6thlc": 0.50,
    "C8": 5.2, "C11": 0.70, "C11b": 0.20, "C11c": 0.85, "C14": 2.5,
    "gamma_coef": 0.12, "gamma_coefb": 0.28, "gamma_coefc": 1.2,
    "mu": 0.0005, "c_K10": 0.35, "c_K10h": 0.35,
    "wpxp_L_thresh": 100.0,
}


def eamv3_params():
    nparams = d.drv_get_nparams()
    params = np.array(d.drv_default_params(nparams))
    idx = {n: i - 1 for n, i in zip(IDX_NAMES, d.drv_param_indices())}
    for name, val in EAMV3_OVERRIDES.items():
        params[idx[name]] = val
    return params, idx


def eam_like_grid(nz=73, ztop=41000.0, seed=None):
    """CLUBB-oriented zi/zt heights (index 1 = surface): EAM-like
    hyperbolic-stretched interfaces (dz ~ 25 m at the surface to
    ~2.5 km at the model top), plus the below-surface ghost point."""
    s = np.linspace(0.0, 1.0, nz)
    zi = ztop * (np.expm1(3.2 * s) / np.expm1(3.2))
    zi[0] = 0.0
    if seed is not None:
        rng = np.random.default_rng(seed)
        mid = zi[1:-1]
        gaps = np.diff(zi)
        jitter = rng.uniform(-0.25, 0.25, mid.size) \
            * np.minimum(gaps[:-1], gaps[1:])
        zi[1:-1] = mid + jitter
    zt = np.empty(nz)
    zt[1:] = 0.5 * (zi[1:] + zi[:-1])
    zt[0] = -zt[1]
    return zi, zt


def gen_grid():
    # nz is fixed at 73 across grids: setup_grid (which sets gr%nz and
    # allocates) runs once inside setup_clubb_core_api, exactly as in
    # EAM where pverp never changes; subsequent drv_setup calls take
    # the setup_grid_heights_api refresh path.
    grids = {"a": eam_like_grid(73), "b": eam_like_grid(73, 28000.0),
             "c": eam_like_grid(73, 41000.0, seed=7)}
    params, _ = eamv3_params()
    out = {}
    rng = np.random.default_rng(42)
    for tag, (zi, zt) in grids.items():
        nz = zi.size
        err = d.drv_setup(29, params, zi, zt)
        assert err == 0
        zm_o, zt_o, dzm, dzt, idzm, idzt, wt2m, wm2t = d.drv_grid_arrays(nz)
        out[f"{tag}_zi_in"] = zi
        out[f"{tag}_zt_in"] = zt
        out[f"{tag}_zm"] = zm_o
        out[f"{tag}_zt"] = zt_o
        out[f"{tag}_dzm"] = dzm
        out[f"{tag}_dzt"] = dzt
        out[f"{tag}_invrs_dzm"] = idzm
        out[f"{tag}_invrs_dzt"] = idzt
        out[f"{tag}_w_zt2zm"] = wt2m
        out[f"{tag}_w_zm2zt"] = wm2t
        # operator sweeps: smooth, oscillatory, random, constant fields
        fields = np.stack([
            300.0 + 0.01 * zt, np.sin(zt / 2000.0),
            rng.normal(size=nz), np.full(nz, 5.0),
            np.exp(-zt / 8000.0)])
        fields_m = np.stack([
            1.0 + zm_o / 1e4, np.cos(zm_o / 3000.0),
            rng.normal(size=nz), np.full(nz, -2.0),
            np.exp(-zm_o / 7000.0)])
        z2m = np.empty_like(fields)
        m2z = np.empty_like(fields)
        dzt_o = np.empty_like(fields)
        dzm_o = np.empty_like(fields)
        for i in range(fields.shape[0]):
            z2m[i], m2z[i], dzt_o[i], dzm_o[i] = \
                d.drv_grid_ops(fields[i], fields_m[i])
        out[f"{tag}_fields_zt"] = fields
        out[f"{tag}_fields_zm"] = fields_m
        out[f"{tag}_zt2zm"] = z2m
        out[f"{tag}_zm2zt"] = m2z
        out[f"{tag}_ddzt"] = dzt_o
        out[f"{tag}_ddzm"] = dzm_o
    return out


def gen_tridag():
    rng = np.random.default_rng(1)
    cases = []
    for n, nrhs in [(73, 1), (73, 2), (73, 5), (41, 1), (5, 3), (2, 1)]:
        # CLUBB-like diagonally dominant (implicit diffusion matrix)
        sub = -rng.uniform(0.1, 1.0, n)
        sup = -rng.uniform(0.1, 1.0, n)
        dia = 1.0 + rng.uniform(0.5, 2.0, n) + np.abs(sub) + np.abs(sup)
        rhs = rng.normal(size=(n, nrhs))
        cases.append((sup, dia, sub, rhs))
        # non-dominant, pivoting exercised
        sub = rng.normal(size=n) * 10.0
        sup = rng.normal(size=n)
        dia = rng.normal(size=n) * 0.1
        rhs = rng.normal(size=(n, nrhs))
        cases.append((sup, dia, sub, rhs))
        # mixed magnitudes
        sub = rng.normal(size=n) * 10.0 ** rng.integers(-6, 6, n)
        sup = rng.normal(size=n) * 10.0 ** rng.integers(-6, 6, n)
        dia = rng.normal(size=n) * 10.0 ** rng.integers(-6, 6, n)
        rhs = rng.normal(size=(n, nrhs))
        cases.append((sup, dia, sub, rhs))
    out = {"ncases": np.array(len(cases))}
    for i, (sup, dia, sub, rhs) in enumerate(cases):
        sol, err = d.drv_tridag_solve(sup, dia, sub, rhs)
        assert err == 0, (i, err)
        out[f"{i}_supd"] = sup
        out[f"{i}_diag"] = dia
        out[f"{i}_subd"] = sub
        out[f"{i}_rhs"] = rhs
        out[f"{i}_solution"] = sol
    return out


def gen_sat():
    t = np.concatenate([np.linspace(150.0, 330.0, 181),
                        [273.15, 273.15 - 85.0, 273.15 - 90.0,
                         273.15 - 84.999, 273.15 - 89.999, 200.0]])
    p = np.concatenate([np.geomspace(50.0, 1.05e5, 181),
                        [1.0e5, 500.0, 120.0, 3000.0, 2.0, 700.0]])
    # all (p, T) pairs, flattened
    pp, tt = np.meshgrid(p, t, indexing="ij")
    pp, tt = pp.ravel(), tt.ravel()
    rsl, rsi = d.drv_sat(pp, tt)
    return {"p": pp, "t": tt, "rsl": rsl, "rsi": rsi}


def _build_columns(nz, zt, params, idx, rng):
    """40 regime-sweeping pdf_closure input columns (CLUBB
    orientation, level 1 = below-surface ghost)."""
    gamma_coef = params[idx["gamma_coef"]]
    gamma_coefb = params[idx["gamma_coefb"]]
    gamma_coefc = params[idx["gamma_coefc"]]
    skw_denom_coef = params[idx["Skw_denom_coef"]]

    ncol = 40
    z = zt.copy()
    cols = {k: np.zeros((ncol, nz)) for k in
            ["p", "exner", "thv_ds", "wm", "wp2", "wp3", "skw", "skthl",
             "skrt", "rtm", "rtp2", "wprtp", "thlm", "thlp2", "wpthlp",
             "um", "up2", "upwp", "vm", "vp2", "vpwp", "rtpthlp",
             "gamma_skw_fnc", "sigma_sqd_w", "wp2_in3", "rtp3", "thlp3"]}

    for c in range(ncol):
        fam = c % 8
        h = rng.uniform(7000.0, 8500.0)
        p = 1.0e5 * np.exp(-np.maximum(z, 0.0) / h)
        p[0] = p[1]
        exner = (p / 1.0e5) ** (287.042 / 1004.64)
        t_sfc = rng.uniform(275.0, 302.0)
        lapse = rng.uniform(0.005, 0.0075)
        T = np.maximum(t_sfc - lapse * np.maximum(z, 0.0), 195.0)
        thlm = T / exner
        rsl, _ = d.drv_sat(p, T)

        zpbl = rng.uniform(800.0, 2500.0)
        bl = np.exp(-np.maximum(z, 0.0) / zpbl)

        if fam == 0:    # stable BL, negative skewness
            rh = rng.uniform(0.3, 0.9)
            wp2 = 0.02 * bl + 1e-3
            skw_t = -rng.uniform(0.2, 1.5) * bl
        elif fam == 1:  # convective, cloudy layer
            rh = rng.uniform(0.85, 1.0)
            wp2 = rng.uniform(0.3, 1.0) * bl + 1e-3
            skw_t = rng.uniform(0.5, 2.5) * bl
        elif fam == 2:  # extreme skewness, both signs
            rh = rng.uniform(0.5, 1.0)
            wp2 = rng.uniform(0.2, 0.8) * bl + 1e-3
            skw_t = np.sign(rng.normal()) * rng.uniform(3.5, 4.4) * bl
        elif fam == 3:  # degenerate wp2 (below w_tol_sqd = 4e-4)
            rh = rng.uniform(0.4, 1.0)
            wp2 = np.full(nz, rng.uniform(1e-6, 3e-4))
            skw_t = np.zeros(nz)
        elif fam == 4:  # degenerate scalar variances in layers
            rh = rng.uniform(0.5, 0.98)
            wp2 = 0.3 * bl + 1e-3
            skw_t = rng.uniform(0.2, 1.0) * bl
        elif fam == 5:  # saturated column (cf ~ 1)
            rh = rng.uniform(1.01, 1.08)
            wp2 = 0.4 * bl + 1e-3
            skw_t = rng.uniform(0.3, 1.5) * bl
        elif fam == 6:  # cold/cirrus emphasis (ice supersat branch)
            rh = rng.uniform(0.7, 1.05)
            wp2 = 0.05 + 0.2 * np.exp(-((z - 9000.0) / 2500.0) ** 2)
            skw_t = 0.8 * np.exp(-((z - 9000.0) / 2500.0) ** 2)
        else:           # randomized
            rh = rng.uniform(0.2, 1.05)
            wp2 = rng.uniform(1e-4, 1.2, nz) * bl + 5e-4
            skw_t = rng.normal(0.0, 1.5, nz) * bl

        rtm = np.clip(rh * rsl, 1e-7, 0.025)
        # moments
        rtp2 = (rng.uniform(0.02, 0.25) * rtm) ** 2 * bl + 1e-14
        thlp2 = rng.uniform(0.05, 1.5) ** 2 * bl + 1e-6
        if fam == 4:
            rtp2[nz // 2:] = 1e-17          # below rt_tol^2
            thlp2[: nz // 4] = 1e-5         # below thl_tol^2
        corr_wrt = rng.uniform(-0.7, 0.9, nz)
        corr_wthl = rng.uniform(-0.9, 0.7, nz)
        corr_rtthl = rng.uniform(-0.95, 0.95, nz)
        wprtp = corr_wrt * np.sqrt(wp2 * rtp2)
        wpthlp = corr_wthl * np.sqrt(wp2 * thlp2)
        rtpthlp = corr_rtthl * np.sqrt(rtp2 * thlp2)
        wp3 = skw_t * wp2 ** 1.5
        rtp3 = rng.normal(0.0, 0.5, nz) * rtp2 ** 1.5
        thlp3 = rng.normal(0.0, 0.5, nz) * thlp2 ** 1.5

        um = 5.0 + 10.0 * np.tanh(z / 3000.0) + rng.normal(0, 2)
        vm = -3.0 + rng.normal(0, 2) * bl
        up2 = rng.uniform(0.05, 0.6) * bl + 1e-4
        vp2 = rng.uniform(0.05, 0.6) * bl + 1e-4
        upwp = rng.uniform(-0.6, 0.6, nz) * np.sqrt(wp2 * up2)
        vpwp = rng.uniform(-0.6, 0.6, nz) * np.sqrt(wp2 * vp2)
        wm = rng.normal(0.0, 0.02, nz) * bl

        thv_ds = thlm.copy()
        # Skewness/sigma exactly as pdf_closure_driver builds them
        skw = d.drv_skx(wp2, wp3, 2.0e-2)
        skthl = d.drv_skx(thlp2, thlp3, 1.0e-2)
        skrt = d.drv_skx(rtp2, rtp3, 1.0e-8)
        gam = gamma_coefb + (gamma_coef - gamma_coefb) * np.exp(
            -0.5 * (skw / gamma_coefc) ** 2)
        sig = d.drv_sigma_sqd_w(gam, wp2, thlp2, rtp2, up2, vp2,
                                wpthlp, wprtp, upwp, vpwp)

        for k, v in [("p", p), ("exner", exner), ("thv_ds", thv_ds),
                     ("wm", wm), ("wp2", wp2), ("wp3", wp3),
                     ("skw", skw), ("skthl", skthl), ("skrt", skrt),
                     ("rtm", rtm), ("rtp2", rtp2), ("wprtp", wprtp),
                     ("thlm", thlm), ("thlp2", thlp2),
                     ("wpthlp", wpthlp), ("um", um), ("up2", up2),
                     ("upwp", upwp), ("vm", vm), ("vp2", vp2),
                     ("vpwp", vpwp), ("rtpthlp", rtpthlp),
                     ("gamma_skw_fnc", gam), ("sigma_sqd_w", sig),
                     ("wp2_in3", wp2), ("rtp3", rtp3),
                     ("thlp3", thlp3)]:
            cols[k][c] = v
    del cols["wp2_in3"]
    assert skw_denom_coef == 0.0  # CLUBB_CAM compile-time default
    return cols


def gen_pdf_closure(params, idx):
    zi, zt = eam_like_grid(73)
    nz = zi.size
    err = d.drv_setup(29, params, zi, zt)
    assert err == 0
    rng = np.random.default_rng(2024)
    cols = _build_columns(nz, zt, params, idx, rng)
    ncol = cols["p"].shape[0]

    moments = np.zeros((ncol, nz, 17))
    pdfp = np.zeros((ncol, nz, 47))
    sig_out = np.zeros((ncol, nz))
    for c in range(ncol):
        m, pp, so, perr = d.drv_pdf_closure(
            cols["p"][c], cols["exner"][c], cols["thv_ds"][c],
            cols["wm"][c], cols["wp2"][c], cols["wp3"][c],
            cols["sigma_sqd_w"][c], cols["skw"][c], cols["skthl"][c],
            cols["skrt"][c], cols["rtm"][c], cols["rtp2"][c],
            cols["wprtp"][c], cols["thlm"][c], cols["thlp2"][c],
            cols["wpthlp"][c], cols["um"][c], cols["up2"][c],
            cols["upwp"][c], cols["vm"][c], cols["vp2"][c],
            cols["vpwp"][c], cols["rtpthlp"][c])
        assert perr == 0, (c, perr)
        moments[c] = m
        pdfp[c] = pp
        sig_out[c] = so
    out = {f"in_{k}": v for k, v in cols.items()}
    out["zi"] = zi
    out["zt"] = zt
    out["moments"] = moments
    out["pdf_params"] = pdfp
    out["sigma_sqd_w_out"] = sig_out
    return out


OUTS_SLOTS = ["rcm", "cloud_frac", "ice_supersat_frac", "wprcp",
              "sigma_sqd_w", "wpthvp", "wp2thvp", "rtpthvp", "thlpthvp",
              "rc_coef", "rcm_in_layer", "cloud_cover", "rcp2_zt",
              "thlprcp", "rc_coef_zm", "wp2rtp", "wp2thlp", "wp2rcp",
              "rtprcp", "rcp2", "uprcp", "vprcp", "cloud_frac_zm",
              "ice_supersat_frac_zm", "rtm_zm", "thlm_zm", "rcm_zm",
              "rcm_supersat_adj", "sigma_sqd_w_zt"]

DRIVER_IN_ZT = ["thlm", "rtm", "rtp3", "thlp3", "wp3", "wm_zt", "um",
                "vm", "p", "exner", "thv_ds_zt", "rfrzm"]
DRIVER_IN_ZM = ["wprtp", "wpthlp", "rtp2", "thlp2", "rtpthlp", "wp2",
                "up2", "upwp", "vp2", "vpwp", "wm_zm", "thv_ds_zm"]
DT_EAMV3 = 300.0  # clubb_timestep


def _interp_zm(z_m, z_t, field_t):
    """Simple linear resample of a zt-built profile onto zm heights
    (input construction only -- NOT CLUBB's zt2zm)."""
    return np.interp(z_m, z_t, field_t)


def _build_driver_columns(nz, zt, zm, rng):
    """40 regime-sweeping pdf_closure_driver input columns with each
    field on its native grid (means/third-moments on zt, second
    moments/fluxes on zm; CLUBB orientation, level 1 = surface/ghost).
    Family 4 includes a sharp dry notch in rtm that drives rcm > rtm
    at the notch after trapezoidal averaging (clip_rcm coverage) and
    cloud-top/base transitions for compute_cloud_cover."""
    ncol = 40
    z_t = np.maximum(zt, 0.0)
    z_m = np.maximum(zm, 0.0)
    cols = {k: np.zeros((ncol, nz)) for k in DRIVER_IN_ZT + DRIVER_IN_ZM}

    for c in range(ncol):
        fam = c % 8
        h = rng.uniform(7000.0, 8500.0)
        p = 1.0e5 * np.exp(-z_t / h)
        p[0] = p[1]
        exner = (p / 1.0e5) ** (287.042 / 1004.64)
        t_sfc = rng.uniform(275.0, 302.0)
        lapse = rng.uniform(0.005, 0.0075)
        T = np.maximum(t_sfc - lapse * z_t, 195.0)
        thlm = T / exner
        rsl, _ = d.drv_sat(p, T)

        zpbl = rng.uniform(800.0, 2500.0)
        bl_t = np.exp(-z_t / zpbl)
        bl_m = np.exp(-z_m / zpbl)

        if fam == 0:    # stable BL, negative skewness
            rh = rng.uniform(0.3, 0.9)
            wp2 = 0.02 * bl_m + 1e-3
            skw_t = -rng.uniform(0.2, 1.5) * bl_t
        elif fam == 1:  # convective, cloudy layer
            rh = rng.uniform(0.85, 1.0)
            wp2 = rng.uniform(0.3, 1.0) * bl_m + 1e-3
            skw_t = rng.uniform(0.5, 2.5) * bl_t
        elif fam == 2:  # extreme skewness, both signs
            rh = rng.uniform(0.5, 1.0)
            wp2 = rng.uniform(0.2, 0.8) * bl_m + 1e-3
            skw_t = np.sign(rng.normal()) * rng.uniform(3.5, 4.4) * bl_t
        elif fam == 3:  # degenerate wp2 (below w_tol_sqd)
            rh = rng.uniform(0.4, 1.0)
            wp2 = np.full(nz, rng.uniform(1e-6, 3e-4))
            skw_t = np.zeros(nz)
        elif fam == 4:  # saturated with a sharp one-level dry notch +
                        # large rt variance: drives the trapezoidal rcm
                        # above the notch rtm (clip_rcm coverage) and
                        # cloud-top/base compute_cloud_cover branches
            rh = rng.uniform(1.0, 1.06)
            wp2 = 0.3 * bl_m + 1e-3
            skw_t = rng.uniform(0.2, 1.0) * bl_t
        elif fam == 5:  # saturated column (cf ~ 1)
            rh = rng.uniform(1.01, 1.08)
            wp2 = 0.4 * bl_m + 1e-3
            skw_t = rng.uniform(0.3, 1.5) * bl_t
        elif fam == 6:  # cold/cirrus emphasis (ice supersat branch)
            rh = rng.uniform(0.7, 1.05)
            wp2 = 0.05 + 0.2 * np.exp(-((z_m - 9000.0) / 2500.0) ** 2)
            skw_t = 0.8 * np.exp(-((z_t - 9000.0) / 2500.0) ** 2)
        else:           # randomized
            rh = rng.uniform(0.2, 1.05)
            wp2 = rng.uniform(1e-4, 1.2, nz) * bl_m + 5e-4
            skw_t = rng.normal(0.0, 1.5, nz) * bl_t

        rtm = np.clip(rh * rsl, 1e-7, 0.025)
        if fam == 4:
            k0 = rng.integers(6, 16)
            rtm[k0] *= rng.uniform(5e-4, 2e-3)  # one-level dry notch

        rtm_m = _interp_zm(z_m, z_t, rtm)
        if fam == 4:
            rtp2 = (rng.uniform(0.4, 0.55) * rtm_m) ** 2 + 1e-14
        else:
            rtp2 = (rng.uniform(0.02, 0.25) * rtm_m) ** 2 * bl_m + 1e-14
        thlp2 = rng.uniform(0.05, 1.5) ** 2 * bl_m + 1e-6
        if fam == 3:
            rtp2[nz // 2:] = 1e-17          # below rt_tol^2
            thlp2[: nz // 4] = 1e-5
        corr_wrt = rng.uniform(-0.7, 0.9, nz)
        corr_wthl = rng.uniform(-0.9, 0.7, nz)
        corr_rtthl = rng.uniform(-0.95, 0.95, nz)
        wprtp = corr_wrt * np.sqrt(wp2 * rtp2)
        wpthlp = corr_wthl * np.sqrt(wp2 * thlp2)
        rtpthlp = corr_rtthl * np.sqrt(rtp2 * thlp2)

        wp2_t = _interp_zm(z_t, z_m, wp2)
        wp3 = skw_t * wp2_t ** 1.5
        rtp2_t = _interp_zm(z_t, z_m, rtp2)
        thlp2_t = _interp_zm(z_t, z_m, thlp2)
        rtp3 = rng.normal(0.0, 0.5, nz) * rtp2_t ** 1.5
        thlp3 = rng.normal(0.0, 0.5, nz) * thlp2_t ** 1.5

        um = 5.0 + 10.0 * np.tanh(z_t / 3000.0) + rng.normal(0, 2)
        vm = -3.0 + rng.normal(0, 2) * bl_t
        up2 = rng.uniform(0.05, 0.6) * bl_m + 1e-4
        vp2 = rng.uniform(0.05, 0.6) * bl_m + 1e-4
        upwp = rng.uniform(-0.6, 0.6, nz) * np.sqrt(wp2 * up2)
        vpwp = rng.uniform(-0.6, 0.6, nz) * np.sqrt(wp2 * vp2)
        wm_zt = rng.normal(0.0, 0.02, nz) * bl_t
        wm_zm = _interp_zm(z_m, z_t, wm_zt)

        thv_ds_zt = thlm.copy()
        thv_ds_zm = _interp_zm(z_m, z_t, thlm)
        rfrzm = np.zeros(nz)

        for k, v in [("thlm", thlm), ("rtm", rtm), ("rtp3", rtp3),
                     ("thlp3", thlp3), ("wp3", wp3), ("wm_zt", wm_zt),
                     ("um", um), ("vm", vm), ("p", p), ("exner", exner),
                     ("thv_ds_zt", thv_ds_zt), ("rfrzm", rfrzm),
                     ("wprtp", wprtp), ("wpthlp", wpthlp),
                     ("rtp2", rtp2), ("thlp2", thlp2),
                     ("rtpthlp", rtpthlp), ("wp2", wp2), ("up2", up2),
                     ("upwp", upwp), ("vp2", vp2), ("vpwp", vpwp),
                     ("wm_zm", wm_zm), ("thv_ds_zm", thv_ds_zm)]:
            cols[k][c] = v
    return cols


def _run_pdf_driver(cols, c):
    args = [cols[k][c] for k in
            ["wprtp", "thlm", "wpthlp", "rtp2", "rtp3", "thlp2",
             "thlp3", "rtpthlp", "wp2", "wp3", "wm_zm", "wm_zt", "um",
             "up2", "upwp", "vm", "vp2", "vpwp", "p", "exner",
             "thv_ds_zm", "thv_ds_zt", "rfrzm", "rtm"]]
    return d.drv_pdf_closure_driver(DT_EAMV3, *args)


def gen_pdf_driver(params, idx):
    zi, zt = eam_like_grid(73)
    nz = zi.size
    err = d.drv_setup(29, params, zi, zt)
    assert err == 0
    d.drv_set_eam_flags(2, True)  # EAMv3: ipdf placement 2, expldiff T

    rng = np.random.default_rng(20260712)
    cols = _build_driver_columns(nz, zt, zi, rng)
    ncol = cols["p"].shape[0]

    rtm_out = np.zeros((ncol, nz))
    outs = np.zeros((ncol, nz, len(OUTS_SLOTS)))
    pdfp_zt = np.zeros((ncol, nz, 47))
    pdfp_zm = np.zeros((ncol, nz, 47))
    for c in range(ncol):
        ro, o, pz, pm, perr = _run_pdf_driver(cols, c)
        assert perr == 0, (c, perr)
        rtm_out[c] = ro
        outs[c] = o
        pdfp_zt[c] = pz
        pdfp_zm[c] = pm
        # l_rtm_nudge = .false.: rtm must pass through unchanged
        assert np.array_equal(ro, cols["rtm"][c]), c

    # Branch coverage sanity: cloud boundaries (compute_cloud_cover's
    # partial-fill branch needs rcm crossing rc_tol) must occur.
    rcm = outs[:, :, OUTS_SLOTS.index("rcm")]
    crossings = ((rcm[:, 1:-1] >= 1e-6)
                 & ((rcm[:, 2:] < 1e-6) | (rcm[:, :-2] < 1e-6))).sum()
    assert crossings > 20, crossings

    # ---- END-TO-END VALIDATION of the extraction + 8 "adv" cases ----
    # Advance the real (public) advance_clubb_core one EAMv3 step from
    # 8 of the synthetic columns; the passthrough outputs of its final
    # (ipdf_call_placement=2) internal pdf_closure_driver call must be
    # reproduced BITWISE by the extracted routine replayed on the
    # advanced state.  The advanced states + outputs are recorded as
    # extra golden cases (realistic covariance structure).
    z_t = np.maximum(zt, 0.0)
    adv_sel = list(range(8))
    adv_in = {k: np.zeros((len(adv_sel), nz)) for k in
              DRIVER_IN_ZT + DRIVER_IN_ZM}
    adv_rtm_out = np.zeros((len(adv_sel), nz))
    adv_outs = np.zeros((len(adv_sel), nz, len(OUTS_SLOTS)))
    adv_pdfp_zt = np.zeros((len(adv_sel), nz, 47))
    adv_pdfp_zm = np.zeros((len(adv_sel), nz, 47))
    ned = 29
    zero = np.zeros(nz)
    passthrough = [("rcm", "prog", 17), ("cloud_frac", "prog", 18),
                   ("wpthvp", "prog", 19), ("wp2thvp", "prog", 20),
                   ("rtpthvp", "prog", 21), ("thlpthvp", "prog", 22),
                   ("rcp2_zt", "diag", 2), ("thlprcp", "diag", 3),
                   ("wprcp", "diag", 4), ("ice_supersat_frac", "diag", 5),
                   ("rcm_in_layer", "diag", 6), ("cloud_cover", "diag", 7)]
    for j, c in enumerate(adv_sel):
        p = cols["p"][c]
        exner = cols["exner"][c]
        T = cols["thlm"][c] * exner
        rho_t = p / (287.042 * T)
        rho_ds_zt = rho_t.copy()
        rho_ds_zt[0] = rho_ds_zt[1]
        invrs_rho_ds_zt = 1.0 / rho_ds_zt
        rho_ds_zm = _interp_zm(np.maximum(zi, 0.0), z_t, rho_ds_zt)
        invrs_rho_ds_zm = 1.0 / rho_ds_zm
        prog_in = np.column_stack(
            [cols["um"][c], cols["vm"][c], cols["upwp"][c],
             cols["vpwp"][c], cols["up2"][c], cols["vp2"][c],
             cols["thlm"][c], cols["rtm"][c], cols["wprtp"][c],
             cols["wpthlp"][c], cols["wp2"][c], cols["wp3"][c],
             cols["rtp2"][c], cols["rtp3"][c], cols["thlp2"][c],
             cols["thlp3"][c], cols["rtpthlp"][c], np.zeros(nz),
             np.zeros(nz), np.zeros(nz), np.zeros(nz), np.zeros(nz),
             np.zeros(nz)])
        edsclr_in = np.zeros((nz, ned))
        for jj in range(ned):
            edsclr_in[:, jj] = (1.0 + 0.1 * jj) * np.exp(
                -z_t / (2000.0 + 300.0 * jj))
        edsclr_in[0, :] = edsclr_in[1, :]
        prog_out, eds_out, diag, apz, apm, aerr = d.drv_advance_clubb_core(
            DT_EAMV3, 1.0e-4, 0.0, 0.02, 5e-5, -0.05, 0.02,
            100000.0, 100000.0,
            zero, zero, zero, zero, zero, zero, zero, zero, zero,
            cols["wm_zm"][c], cols["wm_zt"][c], p, rho_ds_zm, rho_t,
            exner, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm,
            invrs_rho_ds_zt, cols["thv_ds_zm"][c], cols["thv_ds_zt"][c],
            cols["rfrzm"][c], zero, prog_in, edsclr_in)
        assert aerr == 0, (c, aerr)

        # record the advanced state as pdf-driver inputs
        st = {"um": prog_out[:, 0], "vm": prog_out[:, 1],
              "upwp": prog_out[:, 2], "vpwp": prog_out[:, 3],
              "up2": prog_out[:, 4], "vp2": prog_out[:, 5],
              "thlm": prog_out[:, 6], "rtm": prog_out[:, 7],
              "wprtp": prog_out[:, 8], "wpthlp": prog_out[:, 9],
              "wp2": prog_out[:, 10], "wp3": prog_out[:, 11],
              "rtp2": prog_out[:, 12], "rtp3": prog_out[:, 13],
              "thlp2": prog_out[:, 14], "thlp3": prog_out[:, 15],
              "rtpthlp": prog_out[:, 16]}
        for k in DRIVER_IN_ZT + DRIVER_IN_ZM:
            adv_in[k][j] = st.get(k, cols[k][c])

        ro, o, pz, pm, perr = (lambda a: d.drv_pdf_closure_driver(
            DT_EAMV3, a["wprtp"], a["thlm"], a["wpthlp"], a["rtp2"],
            a["rtp3"], a["thlp2"], a["thlp3"], a["rtpthlp"], a["wp2"],
            a["wp3"], cols["wm_zm"][c], cols["wm_zt"][c], a["um"],
            a["up2"], a["upwp"], a["vm"], a["vp2"], a["vpwp"], p,
            exner, cols["thv_ds_zm"][c], cols["thv_ds_zt"][c],
            cols["rfrzm"][c], a["rtm"]))(st)
        assert perr == 0, (c, perr)

        # BITWISE identity between extraction replay and the real
        # advance_clubb_core passthroughs
        for name, kind, slot in passthrough:
            adv_val = prog_out[:, slot] if kind == "prog" \
                else diag[:, slot]
            got = o[:, OUTS_SLOTS.index(name)]
            assert np.array_equal(adv_val, got), (c, name)
        assert np.array_equal(apz, pz), c
        assert np.array_equal(apm, pm), c

        adv_rtm_out[j] = ro
        adv_outs[j] = o
        adv_pdfp_zt[j] = pz
        adv_pdfp_zm[j] = pm
    print("extraction validated bitwise against advance_clubb_core "
          f"on {len(adv_sel)} columns")

    out = {f"in_{k}": v for k, v in cols.items()}
    out.update({f"adv_in_{k}": v for k, v in adv_in.items()})
    out.update(zi=zi, zt=zt, dt=np.array(DT_EAMV3),
               rtm_out=rtm_out, outs=outs, pdfp_zt=pdfp_zt,
               pdfp_zm=pdfp_zm, adv_rtm_out=adv_rtm_out,
               adv_outs=adv_outs, adv_pdfp_zt=adv_pdfp_zt,
               adv_pdfp_zm=adv_pdfp_zm)
    return out


# ---------------------------------------------------------------------------
# Slice C: advance_xp2_xpyp + clip_covars_denom
# ---------------------------------------------------------------------------

XP2_IDX_NAMES = ["C4", "C5", "C14", "c_K2", "c_K9", "nu2", "nu9"]

# Everything advance_xp2_xpyp takes, on its native grid (zm = momentum,
# zt = thermodynamic).  wprtp2/wpthlp2/wprtpthlp are real arguments but
# UNUSED under EAMv3 (l_explicit_turbulent_adv_xpyp=F); recorded as
# zeros.  skw_zm and cloud_frac are likewise dead (l_single_C2_Skw=F,
# l_C2_cloud_frac=F) but recorded with consistent values.
XP2_INPUTS = ["tau_zm", "wm_zm", "rtm", "wprtp", "thlm", "wpthlp",
              "wpthvp", "um", "vm", "wp2", "wp2_zt", "wp3", "upwp",
              "vpwp", "sigma_sqd_w", "skw_zm", "wprtp2", "wpthlp2",
              "wprtpthlp", "kh_zt", "rtp2_forcing", "thlp2_forcing",
              "rtpthlp_forcing", "rho_ds_zm", "rho_ds_zt",
              "invrs_rho_ds_zm", "thv_ds_zm", "cloud_frac", "lscale",
              "wp3_on_wp2", "wp3_on_wp2_zt", "wp2_splat",
              "rtp2", "thlp2", "rtpthlp", "up2", "vp2"]

XP2_PROG = ["rtp2", "thlp2", "rtpthlp", "up2", "vp2"]

W_TOL_SQD = 4.0e-4
RT_TOL = 1.0e-8
THL_TOL = 1.0e-2


def _zt2zm_f(azt, nz):
    return d.drv_grid_ops(azt, np.zeros(nz))[0]


def _zm2zt_f(azm, nz):
    return d.drv_grid_ops(np.zeros(nz), azm)[1]


def _xp2_env(nz, z_t, z_m, rng):
    """Plausible thermodynamic environment + tau/Kh/Lscale profiles
    (inputs to advance_xp2_xpyp that its validation replays verbatim;
    physical consistency only sets regime coverage)."""
    h = rng.uniform(7000.0, 8500.0)
    p = 1.0e5 * np.exp(-z_t / h)
    p[0] = p[1]
    exner = (p / 1.0e5) ** (287.042 / 1004.64)
    t_sfc = rng.uniform(275.0, 302.0)
    T = np.maximum(t_sfc - rng.uniform(0.005, 0.0075) * z_t, 195.0)
    thlm = T / exner
    rsl, _ = d.drv_sat(p, T)
    rtm = np.clip(rng.uniform(0.3, 1.05) * rsl, 1e-7, 0.025)
    rho_ds_zt = p / (287.042 * T)
    rho_ds_zt[0] = rho_ds_zt[1]
    rho_ds_zm = _interp_zm(z_m, z_t, rho_ds_zt)
    zpbl = rng.uniform(800.0, 2500.0)
    bl_m = np.exp(-z_m / zpbl)
    bl_t = np.exp(-z_t / zpbl)
    tau_zm = rng.uniform(50.0, 400.0) + rng.uniform(400.0, 3000.0) * (
        1.0 - np.exp(-z_m / 1500.0))
    kh_zt = rng.uniform(1.0, 80.0) * bl_t + 0.1
    lscale = rng.uniform(30.0, 1200.0) * bl_t + 10.0
    um = 5.0 + 10.0 * np.tanh(z_t / 3000.0) + rng.normal(0, 2)
    vm = -3.0 + rng.normal(0, 2) * bl_t
    wm_zm = rng.normal(0.0, 0.02) * bl_m
    return dict(p=p, thlm=thlm, rtm=rtm, rho_ds_zt=rho_ds_zt,
                rho_ds_zm=rho_ds_zm, bl_m=bl_m, bl_t=bl_t,
                tau_zm=tau_zm, kh_zt=kh_zt, lscale=lscale, um=um,
                vm=vm, wm_zm=wm_zm)


def _finish_xp2_case(nz, env, wp2, wp3, moments, rng):
    """Assemble the full advance_xp2_xpyp input dict from an
    environment, wp2/wp3, and the second-moment dict."""
    wp2 = np.maximum(wp2, W_TOL_SQD)
    wp2_zt = np.maximum(_zm2zt_f(wp2, nz), W_TOL_SQD)
    wp3_zm = _zt2zm_f(wp3, nz)
    wp3_on_wp2_zt = wp3 / wp2_zt
    wp3_on_wp2 = wp3_zm / wp2
    skw_zm = d.drv_skx(wp2, wp3_zm, 2.0e-2)
    gam = 0.28 + (0.12 - 0.28) * np.exp(-0.5 * (skw_zm / 1.2) ** 2)
    sig = d.drv_sigma_sqd_w(
        gam, wp2, moments["thlp2"], moments["rtp2"], moments["up2"],
        moments["vp2"], moments["wpthlp"], moments["wprtp"],
        moments["upwp"], moments["vpwp"])
    case = dict(
        tau_zm=env["tau_zm"], wm_zm=env["wm_zm"], rtm=env["rtm"],
        thlm=env["thlm"], um=env["um"], vm=env["vm"],
        wp2=wp2, wp2_zt=wp2_zt, wp3=wp3, sigma_sqd_w=sig,
        skw_zm=skw_zm, wprtp2=np.zeros(nz), wpthlp2=np.zeros(nz),
        wprtpthlp=np.zeros(nz), kh_zt=env["kh_zt"],
        rtp2_forcing=np.zeros(nz), thlp2_forcing=np.zeros(nz),
        rtpthlp_forcing=np.zeros(nz), rho_ds_zm=env["rho_ds_zm"],
        rho_ds_zt=env["rho_ds_zt"],
        invrs_rho_ds_zm=1.0 / env["rho_ds_zm"],
        thv_ds_zm=_interp_zm(np.maximum(env.get("z_m"), 0.0),
                             np.maximum(env.get("z_t"), 0.0),
                             env["thlm"]),
        cloud_frac=env.get("cloud_frac", np.zeros(nz)),
        lscale=env["lscale"], wp3_on_wp2=wp3_on_wp2,
        wp3_on_wp2_zt=wp3_on_wp2_zt, wp2_splat=np.zeros(nz),
        **{k: moments[k] for k in ["wprtp", "wpthlp", "wpthvp",
                                   "upwp", "vpwp", "rtp2", "thlp2",
                                   "rtpthlp", "up2", "vp2"]})
    del rng
    return case


def _build_xp2_synthetic(nz, zt, zm, rng, ncase=24):
    """Synthetic stress columns for advance_xp2_xpyp."""
    z_t = np.maximum(zt, 0.0)
    z_m = np.maximum(zm, 0.0)
    cases = []
    for c in range(ncase):
        fam = c % 8
        env = _xp2_env(nz, z_t, z_m, rng)
        env["z_t"], env["z_m"] = z_t, z_m
        bl_m, bl_t = env["bl_m"], env["bl_t"]
        rtm_m = _interp_zm(z_m, z_t, env["rtm"])

        wp2 = rng.uniform(0.2, 0.8) * bl_m + 1e-3
        skw_t = rng.normal(0.0, 1.0, nz) * bl_t
        rtp2 = (rng.uniform(0.05, 0.25) * rtm_m) ** 2 * bl_m + 1e-14
        thlp2 = rng.uniform(0.1, 1.0) ** 2 * bl_m + 1e-6
        up2 = rng.uniform(0.05, 0.6) * bl_m + 1e-4
        vp2 = rng.uniform(0.05, 0.6) * bl_m + 1e-4
        corr = lambda lo, hi: rng.uniform(lo, hi, nz)  # noqa: E731

        if fam == 0:
            # degenerate variances: every floor clip + dp1 threshold
            wp2 = np.full(nz, W_TOL_SQD)
            rtp2 = np.full(nz, 1e-18)
            thlp2 = np.full(nz, 1e-9)
            up2 = np.full(nz, 1e-6)
            vp2 = np.full(nz, 1e-6)
            skw_t = np.zeros(nz)
        elif fam == 1:
            # Cauchy-Schwarz-violating rtpthlp (clip_covar) and
            # correlation-violating fluxes (clip_covars_denom bait)
            pass  # handled below via corr multipliers
        elif fam == 2:
            # strong skewness, sign-alternating wp3 (both upwind arms)
            skw_t = 3.5 * np.sign(np.sin(z_t / 900.0)) * bl_t \
                + rng.normal(0.0, 0.3, nz) * bl_t
        elif fam == 3:
            # large rtp2 above the 0.5*rtm^2 cap + positive forcing
            rtp2 = (0.9 * rtm_m) ** 2 + 1e-14
        elif fam == 4:
            # near/above the 1000 m2/s2 up2/vp2 cap + strong shear
            up2 = 900.0 + 600.0 * bl_m
            vp2 = 800.0 + 500.0 * bl_m
            wp2 = 1.5 * bl_m + 1e-3
            env["um"] = env["um"] + 25.0 * np.tanh((z_t - 1500.0) / 300.0)
        elif fam == 5:
            # hole bait: strongly negative forcings drive the solution
            # negative in a band -> pos_definite_variances fires
            pass  # forcings set below
        elif fam == 6:
            # advection/dissipation heavy
            env["wm_zm"] = 0.2 * np.sin(z_m / 2500.0)
            env["tau_zm"] = np.full(nz, rng.uniform(20.0, 60.0))
            env["kh_zt"] = np.full(nz, rng.uniform(80.0, 200.0))
        # fam 7: randomized defaults above

        wp3 = skw_t * np.maximum(_zm2zt_f(wp2, nz), W_TOL_SQD) ** 1.5

        cw = 2.5 if fam == 1 else 0.9
        moments = dict(
            rtp2=rtp2, thlp2=thlp2, up2=up2, vp2=vp2,
            rtpthlp=(corr(-cw, cw) if fam != 1 else
                     np.where(np.arange(nz) % 2 == 0, 2.5, -2.5))
            * np.sqrt(rtp2 * thlp2),
            wprtp=corr(-0.7, 0.9) * np.sqrt(wp2 * rtp2)
            * (1.6 if fam == 1 else 1.0),
            wpthlp=corr(-0.9, 0.7) * np.sqrt(wp2 * thlp2)
            * (1.6 if fam == 1 else 1.0),
            upwp=corr(-0.6, 0.6) * np.sqrt(wp2 * up2)
            * (1.8 if fam == 1 else 1.0),
            vpwp=corr(-0.6, 0.6) * np.sqrt(wp2 * vp2)
            * (1.8 if fam == 1 else 1.0),
            wpthvp=rng.normal(0.0, 0.05, nz) * bl_m)

        case = _finish_xp2_case(nz, env, wp2, wp3, moments, rng)

        if fam == 3:
            case["rtp2_forcing"] = 5e-9 * bl_m
        if fam == 5:
            band = np.exp(-((z_m - 2000.0) / 500.0) ** 2)
            case["rtp2_forcing"] = -3.0 * case["rtp2"] / DT_EAMV3 * band
            case["thlp2_forcing"] = -3.0 * case["thlp2"] / DT_EAMV3 * band
            case["rtpthlp_forcing"] = 2.0 * case["rtpthlp"] / DT_EAMV3 \
                * band
        if fam == 4:
            case["wp2_splat"] = -0.02 * case["wp2"] / DT_EAMV3
        cases.append(case)
    return cases


def _build_xp2_adv(nz, zt, zm, rng, nadv=16):
    """Realistic cases: prognostic moments and Kh from REAL one-step
    advance_clubb_core-advanced states of the slice-B regime columns
    (tau/Lscale profiles are plausible synthetics -- their computation
    is slice D; advance_xp2_xpyp validation replays whatever it is
    fed)."""
    z_t = np.maximum(zt, 0.0)
    z_m = np.maximum(zm, 0.0)
    cols = _build_driver_columns(nz, zt, zm, rng)
    ned = 29
    zero = np.zeros(nz)
    cases = []
    for c in range(nadv):
        p = cols["p"][c]
        exner = cols["exner"][c]
        T = cols["thlm"][c] * exner
        rho_t = p / (287.042 * T)
        rho_ds_zt = rho_t.copy()
        rho_ds_zt[0] = rho_ds_zt[1]
        rho_ds_zm = _interp_zm(z_m, z_t, rho_ds_zt)
        prog_in = np.column_stack(
            [cols["um"][c], cols["vm"][c], cols["upwp"][c],
             cols["vpwp"][c], cols["up2"][c], cols["vp2"][c],
             cols["thlm"][c], cols["rtm"][c], cols["wprtp"][c],
             cols["wpthlp"][c], cols["wp2"][c], cols["wp3"][c],
             cols["rtp2"][c], cols["rtp3"][c], cols["thlp2"][c],
             cols["thlp3"][c], cols["rtpthlp"][c]]
            + [np.zeros(nz)] * 6)
        edsclr_in = np.zeros((nz, ned))
        for jj in range(ned):
            edsclr_in[:, jj] = (1.0 + 0.1 * jj) * np.exp(
                -z_t / (2000.0 + 300.0 * jj))
        edsclr_in[0, :] = edsclr_in[1, :]
        prog_out, _eds, diag, _pz, _pm, aerr = d.drv_advance_clubb_core(
            DT_EAMV3, 1.0e-4, 0.0, 0.02, 5e-5, -0.05, 0.02,
            100000.0, 100000.0,
            zero, zero, zero, zero, zero, zero, zero, zero, zero,
            cols["wm_zm"][c], cols["wm_zt"][c], p, rho_ds_zm, rho_t,
            exner, rho_ds_zm, rho_ds_zt, 1.0 / rho_ds_zm,
            1.0 / rho_ds_zt, cols["thv_ds_zm"][c], cols["thv_ds_zt"][c],
            cols["rfrzm"][c], zero, prog_in, edsclr_in)
        assert aerr == 0, (c, aerr)

        env = _xp2_env(nz, z_t, z_m, rng)
        env["z_t"], env["z_m"] = z_t, z_m
        # replace environment pieces with the real column/advanced state
        env["thlm"] = prog_out[:, 6]
        env["rtm"] = prog_out[:, 7]
        env["um"] = prog_out[:, 0]
        env["vm"] = prog_out[:, 1]
        env["rho_ds_zt"] = rho_ds_zt
        env["rho_ds_zm"] = rho_ds_zm
        env["wm_zm"] = cols["wm_zm"][c]
        env["kh_zt"] = np.maximum(diag[:, 1], 0.0) + 0.1  # real khzt
        env["cloud_frac"] = prog_out[:, 18]
        moments = dict(
            rtp2=prog_out[:, 12], thlp2=prog_out[:, 14],
            rtpthlp=prog_out[:, 16], up2=prog_out[:, 4],
            vp2=prog_out[:, 5], wprtp=prog_out[:, 8],
            wpthlp=prog_out[:, 9], upwp=prog_out[:, 2],
            vpwp=prog_out[:, 3], wpthvp=prog_out[:, 19])
        cases.append(_finish_xp2_case(nz, env, prog_out[:, 10],
                                      prog_out[:, 11], moments, rng))
    return cases


def _run_xp2(case):
    args = [case[k] for k in XP2_INPUTS]
    return d.drv_advance_xp2_xpyp(DT_EAMV3, True, *args)


def gen_xp2(params, idx):
    zi, zt = eam_like_grid(73)
    nz = zi.size
    err = d.drv_setup(29, params, zi, zt)
    assert err == 0
    d.drv_set_eam_flags(2, True)

    nu2, nu9, iflags = d.drv_xp2_config(nz)
    # EAMv3 model_flags configuration (see drv_xp2_config for slots):
    # l_single_C2_Skw=F, l_explicit_turbulent_adv_xpyp=F,
    # l_upwind_xpyp_ta=T, l_min_xp2_from_corr_wx=F, l_C2_cloud_frac=F,
    # l_hole_fill=T, l_tke_aniso=T, up2_vp2 sponge damping OFF.
    assert list(iflags) == [0, 0, 1, 0, 0, 1, 1, 0], list(iflags)

    rng = np.random.default_rng(20260713)
    cases = _build_xp2_adv(nz, zt, zi, rng) \
        + _build_xp2_synthetic(nz, zt, zi, rng)
    ncase = len(cases)

    ins = {k: np.zeros((ncase, nz)) for k in XP2_INPUTS}
    outs = {k: np.zeros((ncase, nz)) for k in XP2_PROG}
    for c, case in enumerate(cases):
        for k in XP2_INPUTS:
            ins[k][c] = case[k]
        r2, t2, rt, u2, v2, xerr = _run_xp2(case)
        assert xerr == 0, (c, xerr)
        outs["rtp2"][c] = r2
        outs["thlp2"][c] = t2
        outs["rtpthlp"][c] = rt
        outs["up2"][c] = u2
        outs["vp2"][c] = v2

    # ---- branch-coverage sanity ----
    assert (outs["rtp2"] == RT_TOL ** 2).any()            # floor clip
    assert (outs["thlp2"] == THL_TOL ** 2).any()          # floor clip
    assert (outs["up2"] == 1000.0).any()                  # magnitude cap
    cap = 0.5 * ins["rtm"] ** 2
    assert (outs["rtp2"] == cap).any()                    # l_clip_large_rtp2
    bound = 0.99 * np.sqrt(outs["rtp2"] * outs["thlp2"])
    at_bound = np.isclose(np.abs(outs["rtpthlp"]), bound, rtol=1e-14) \
        & (bound > 0)
    assert at_bound.any()                                 # clip_covar
    print(f"xp2 coverage: rtp2 floor {(outs['rtp2'] == RT_TOL**2).sum()}, "
          f"cap {(outs['rtp2'] == cap).sum()}, "
          f"up2 cap {(outs['up2'] == 1000.0).sum()}, "
          f"rtpthlp clipped {at_bound.sum()}")

    # ---- clip_covars_denom (post-xp2 instance: cl_num 2/2/2/1/1; the
    # cl_num values only gate stats, which are off) ----
    cc_ins = {k: np.zeros((ncase, nz)) for k in
              ["wp2", "rtp2", "thlp2", "up2", "vp2", "wprtp", "wpthlp",
               "upwp", "vpwp"]}
    cc_outs = {k: np.zeros((ncase, nz)) for k in
               ["wprtp", "wpthlp", "upwp", "vpwp"]}
    for c, case in enumerate(cases):
        cc = dict(wp2=case["wp2"], rtp2=outs["rtp2"][c],
                  thlp2=outs["thlp2"][c], up2=outs["up2"][c],
                  vp2=outs["vp2"][c], wprtp=case["wprtp"],
                  wpthlp=case["wpthlp"], upwp=case["upwp"],
                  vpwp=case["vpwp"])
        wr, wt, uw, vw, cerr = d.drv_clip_covars_denom(
            DT_EAMV3, cc["rtp2"], cc["thlp2"], cc["up2"], cc["vp2"],
            cc["wp2"], 2, 2, 2, 1, 1, cc["wprtp"], cc["wpthlp"],
            cc["upwp"], cc["vpwp"])
        assert cerr == 0, (c, cerr)
        for k in cc_ins:
            cc_ins[k][c] = cc[k]
        cc_outs["wprtp"][c] = wr
        cc_outs["wpthlp"][c] = wt
        cc_outs["upwp"][c] = uw
        cc_outs["vpwp"][c] = vw
    nchanged = sum((cc_outs[k] != cc_ins[k]).any(axis=1).sum()
                   for k in cc_outs)
    assert nchanged > 10, nchanged  # clipping must actually fire
    print(f"clip_covars_denom coverage: {nchanged} clipped columns")

    out = {f"in_{k}": v for k, v in ins.items()}
    out.update({f"out_{k}": v for k, v in outs.items()})
    out.update({f"cc_in_{k}": v for k, v in cc_ins.items()})
    out.update({f"cc_out_{k}": v for k, v in cc_outs.items()})
    out.update(zi=zi, zt=zt, dt=np.array(DT_EAMV3),
               nu2_vert_res_dep=nu2, nu9_vert_res_dep=nu9,
               xp2_flags=np.asarray(iflags),
               nadv=np.array(16))
    return out


# ---------------------------------------------------------------------------
# Slice D: Lscale/tau infrastructure
# ---------------------------------------------------------------------------

LSCALE_IDX_NAMES = ["c_K", "taumin", "taumax", "Lscale_mu_coef",
                    "Lscale_pert_coef", "lmin_coef", "C_wp2_splat",
                    "lambda0_stability_coef", "up2_vp2_factor"]

# drv_lscale_tau_segment per-column array inputs, in signature order
SEG_IN = ["thlm", "rtm", "rcm", "wp2", "wp3", "up2", "vp2", "um", "vm",
          "p", "exner", "thv_ds_zt"]
SEG_SLOTS = ["em", "thvm", "sqrt_em_zt", "lscale", "lscale_up",
             "lscale_down", "tau_zt", "tau_zm", "kh_zt", "kh_zm",
             "wp2_splat", "wp3_splat", "stability_correction",
             "tau_n2_zm"]
SFC_SLOTS = ["wp2", "up2", "vp2", "thlp2", "rtp2", "rtpthlp"]

ML_IN = ["thvm", "thlm", "rtm", "em", "p", "exner", "thv_ds"]
BV_IN = ["thlm", "exner", "rtm", "rcm", "p", "thvm"]

HOST_DXY = 100000.0
MU_EAMV3 = 0.0005

# surface-flux sweep for the segment cases (wpthlp, wprtp, upwp, vpwp):
# stable/unstable, moisture-flux signs, near-zero and strong stress
SEG_SFC_FLUXES = [
    (0.02, 5e-5, -0.05, 0.02),     # weakly convective (slice B/C values)
    (-0.05, -2e-5, -0.02, 0.01),   # stable BL, downward moisture
    (0.3, 8e-4, -0.4, 0.3),        # strongly convective, strong stress
    (-0.3, 2e-4, 0.15, -0.2),      # strongly stable, opposing stress
    (1e-13, 1e-14, 1e-8, -1e-8),   # near-zero everything (ufmin floor)
    (0.0, 0.0, 0.0, 0.0),          # exactly zero
    (0.08, -5e-5, -0.6, -0.5),     # convective, drying, strong shear
    (-0.02, 6e-4, 0.0, 0.0),       # zero stress, mixed scalar fluxes
]


def _advance_entry(cols, c, nz, z_t, z_m, fluxes, sfc_elevation=0.0):
    """Run drv_advance_clubb_core one EAMv3 step from synthetic column
    c (the slice-B pattern); returns (entry_state_dict, prog_out,
    diag)."""
    ned = 29
    zero = np.zeros(nz)
    p = cols["p"][c]
    exner = cols["exner"][c]
    T = cols["thlm"][c] * exner
    rho_t = p / (287.042 * T)
    rho_ds_zt = rho_t.copy()
    rho_ds_zt[0] = rho_ds_zt[1]
    rho_ds_zm = _interp_zm(z_m, z_t, rho_ds_zt)
    prog_in = np.column_stack(
        [cols["um"][c], cols["vm"][c], cols["upwp"][c],
         cols["vpwp"][c], cols["up2"][c], cols["vp2"][c],
         cols["thlm"][c], cols["rtm"][c], cols["wprtp"][c],
         cols["wpthlp"][c], cols["wp2"][c], cols["wp3"][c],
         cols["rtp2"][c], cols["rtp3"][c], cols["thlp2"][c],
         cols["thlp3"][c], cols["rtpthlp"][c]] + [np.zeros(nz)] * 6)
    edsclr_in = np.zeros((nz, ned))
    for jj in range(ned):
        edsclr_in[:, jj] = (1.0 + 0.1 * jj) * np.exp(
            -z_t / (2000.0 + 300.0 * jj))
    edsclr_in[0, :] = edsclr_in[1, :]
    wpthlp_s, wprtp_s, upwp_s, vpwp_s = fluxes
    prog_out, _eds, diag, _pz, _pm, aerr = d.drv_advance_clubb_core(
        DT_EAMV3, 1.0e-4, sfc_elevation, wpthlp_s, wprtp_s, upwp_s,
        vpwp_s, HOST_DXY, HOST_DXY,
        zero, zero, zero, zero, zero, zero, zero, zero, zero,
        cols["wm_zm"][c], cols["wm_zt"][c], p, rho_ds_zm, rho_t,
        exner, rho_ds_zm, rho_ds_zt, 1.0 / rho_ds_zm,
        1.0 / rho_ds_zt, cols["thv_ds_zm"][c], cols["thv_ds_zt"][c],
        cols["rfrzm"][c], zero, prog_in, edsclr_in)
    assert aerr == 0, (c, aerr)
    return prog_in, prog_out, diag, dict(
        p=p, exner=exner, rho_ds_zm=rho_ds_zm, rho_ds_zt=rho_ds_zt,
        rho_t=rho_t, edsclr_in=edsclr_in)


def _seg_case_from_state(cols, c, st):
    """Assemble drv_lscale_tau_segment inputs from a prognostic state
    dict (thlm/rtm/rcm/wp2/wp3/up2/vp2/um/vm) + column environment."""
    return {"thlm": st["thlm"], "rtm": st["rtm"], "rcm": st["rcm"],
            "wp2": st["wp2"], "wp3": st["wp3"], "up2": st["up2"],
            "vp2": st["vp2"], "um": st["um"], "vm": st["vm"],
            "p": cols["p"][c], "exner": cols["exner"][c],
            "thv_ds_zt": cols["thv_ds_zt"][c]}


def gen_lscale(params, idx):
    zi, zt = eam_like_grid(73)
    nz = zi.size
    err = d.drv_setup(29, params, zi, zt)
    assert err == 0
    d.drv_set_eam_flags(2, True)
    d.drv_set_bv_flags(False, False)      # EAMv3 defaults

    lmin, t0, iflags = d.drv_lscale_config()
    # EAMv3 model_flags configuration (drv_lscale_config slots):
    # l_stability_correct_tau_zm=T, l_diag_Lscale_from_tau=F,
    # l_use_C7_Richardson=F, l_use_C11_Richardson=F, l_use_wp3_pr3=F
    # (=> Cx_fnc_Richardson = 0, dead), l_Lscale_plume_centered=F,
    # l_use_ice_latent=F, l_brunt_vaisala_freq_moist=F,
    # l_use_thvm_in_bv_freq=F, l_sat_mixrat_lookup=F, l_tke_aniso=T.
    assert list(iflags) == [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], list(iflags)
    assert lmin == 4.0 and t0 == 300.0, (lmin, t0)
    lidx = {n: int(i) - 1 for n, i in
            zip(LSCALE_IDX_NAMES, d.drv_param_indices_lscale())}
    assert params[lidx["C_wp2_splat"]] == 0.0    # EAMv3 (no override)
    assert params[lidx["c_K"]] == 0.2

    rng = np.random.default_rng(20260714)
    cols = _build_driver_columns(nz, zt, zi, rng)
    z_t = np.maximum(zt, 0.0)
    z_m = np.maximum(zi, 0.0)

    out = {"zi": zi, "zt": zt, "dt": np.array(DT_EAMV3),
           "lmin": np.array(lmin), "T0": np.array(t0),
           "lscale_flags": np.asarray(iflags),
           "host_dxy": np.array(HOST_DXY)}

    # ---- SEG: the advance_clubb_core Lscale/tau segment -------------
    # 16 entry states (rcm=0) + 8 one-step-advanced states (real rcm
    # etc.) + 1 elevated-surface case (else branch).  Validated at
    # generation: kh_zt/kh_zm from the transcription must match the
    # khzt/khzm diagnostics of the REAL advance_clubb_core BITWISE on
    # the same entry state.
    nseg = 25
    seg_ins = {k: np.zeros((nseg, nz)) for k in SEG_IN}
    seg_fluxes = np.zeros((nseg, 4))
    seg_elev = np.zeros(nseg)
    seg_outs = np.zeros((nseg, nz, len(SEG_SLOTS)))
    sfc_outs = np.zeros((nseg, len(SFC_SLOTS)))
    for j in range(nseg):
        c = j % 16
        fluxes = SEG_SFC_FLUXES[j % len(SEG_SFC_FLUXES)]
        elev = 250.0 if j == 24 else 0.0
        prog_in, prog_out, diag, env = _advance_entry(
            cols, c, nz, z_t, z_m, fluxes, sfc_elevation=elev)
        if 16 <= j < 24:
            # advanced state (realistic covariances + non-zero rcm)
            st = dict(um=prog_out[:, 0], vm=prog_out[:, 1],
                      up2=prog_out[:, 4], vp2=prog_out[:, 5],
                      thlm=prog_out[:, 6], rtm=prog_out[:, 7],
                      wp2=prog_out[:, 10], wp3=prog_out[:, 11],
                      rcm=prog_out[:, 17])
        else:
            st = dict(um=prog_in[:, 0], vm=prog_in[:, 1],
                      up2=prog_in[:, 4], vp2=prog_in[:, 5],
                      thlm=prog_in[:, 6], rtm=prog_in[:, 7],
                      wp2=prog_in[:, 10], wp3=prog_in[:, 11],
                      rcm=prog_in[:, 17])
        case = _seg_case_from_state(cols, c, st)
        seg, sfc, serr = d.drv_lscale_tau_segment(
            DT_EAMV3, elev, HOST_DXY, HOST_DXY, *fluxes,
            *[case[k] for k in SEG_IN])
        assert serr == 0, (j, serr)

        # END-TO-END VALIDATION of the transcription: run the REAL
        # advance_clubb_core from this exact state; its khzt/khzm
        # diagnostics are assigned straight from the segment's
        # Kh_zt/Kh_zm before any prognostic advance.
        if 16 <= j < 24:
            prog_in2 = prog_in.copy()
            for slot, key in [(0, "um"), (1, "vm"), (4, "up2"),
                              (5, "vp2"), (6, "thlm"), (7, "rtm"),
                              (10, "wp2"), (11, "wp3"), (17, "rcm")]:
                prog_in2[:, slot] = st[key]
            for slot, col in [(2, 2), (3, 3), (8, 8), (9, 9), (12, 12),
                              (13, 13), (14, 14), (15, 15), (16, 16)]:
                prog_in2[:, slot] = prog_out[:, col]
            zero = np.zeros(nz)
            _po, _e, diag2, _z, _m2, aerr = d.drv_advance_clubb_core(
                DT_EAMV3, 1.0e-4, elev, *fluxes, HOST_DXY, HOST_DXY,
                zero, zero, zero, zero, zero, zero, zero, zero, zero,
                cols["wm_zm"][c], cols["wm_zt"][c], env["p"],
                env["rho_ds_zm"], env["rho_t"], env["exner"],
                env["rho_ds_zm"], env["rho_ds_zt"],
                1.0 / env["rho_ds_zm"], 1.0 / env["rho_ds_zt"],
                cols["thv_ds_zm"][c], cols["thv_ds_zt"][c],
                cols["rfrzm"][c], zero, prog_in2, env["edsclr_in"])
            assert aerr == 0, (j, aerr)
            diag = diag2
        assert np.array_equal(seg[:, 9], diag[:, 0]), j   # kh_zm
        assert np.array_equal(seg[:, 8], diag[:, 1]), j   # kh_zt

        for k in SEG_IN:
            seg_ins[k][j] = case[k]
        seg_fluxes[j] = fluxes
        seg_elev[j] = elev
        seg_outs[j] = seg
        sfc_outs[j] = sfc
    print(f"lscale segment transcription validated bitwise (kh_zm/"
          f"kh_zt) against advance_clubb_core on {nseg} cases")
    # the elevated-surface case must hit the else branch (tolerances)
    assert sfc_outs[24, 3] == THL_TOL ** 2 and sfc_outs[24, 5] == 0.0

    out.update({f"seg_in_{k}": v for k, v in seg_ins.items()})
    out.update(seg_fluxes=seg_fluxes, seg_sfc_elevation=seg_elev,
               seg_outs=seg_outs, seg_sfc_outs=sfc_outs)

    # ---- ML: compute_mixing_length direct kernel sweeps -------------
    nml = 40
    ml_ins = {k: np.zeros((nml, nz)) for k in ML_IN}
    ml_mu = np.zeros(nml)
    ml_lmax = np.zeros(nml)
    ml_outs = {k: np.zeros((nml, nz)) for k in
               ["lscale", "lscale_up", "lscale_down"]}
    for c in range(nml):
        fam = c % 8
        mu = MU_EAMV3
        lmax = 0.25 * HOST_DXY
        if fam == 0:
            # physically consistent: the segment's own thvm/em
            j = c // 8 * 3 % 25
            thvm = seg_outs[j, :, 1]
            em = seg_outs[j, :, 0]
            thlm = seg_ins["thlm"][j]
            rtm = seg_ins["rtm"][j]
            p = seg_ins["p"][j]
            exner = seg_ins["exner"][j]
            thv_ds = seg_ins["thv_ds_zt"][j]
        else:
            h = rng.uniform(7000.0, 8500.0)
            p = 1.0e5 * np.exp(-z_t / h)
            p[0] = p[1]
            exner = (p / 1.0e5) ** (287.042 / 1004.64)
            zpbl = rng.uniform(600.0, 2200.0)
            bl_m = np.exp(-z_m / zpbl)
            if fam == 1:
                # strongly stable BL: sharp inversion, weak TKE ->
                # early parcel exhaustion (quadratic branches)
                T = 285.0 + 0.008 * z_t - 0.004 * np.maximum(
                    z_t - 3000.0, 0.0)
                thlm = T / exner
                rtm = np.full(nz, rng.uniform(1e-4, 2e-3))
                em = rng.uniform(0.005, 0.05) * bl_m + 6.0e-4
            elif fam == 2:
                # strongly unstable: superadiabatic thvm, huge TKE ->
                # parcels run far; a small host-grid Lscale_max
                # (host_dx = 8 km) makes the cap branch fire
                T = rng.uniform(295.0, 305.0) - 0.0099 * z_t
                thlm = np.maximum(T, 210.0) / exner
                thlm = np.minimum.accumulate(thlm + 1e-3 * z_t) \
                    - 1e-3 * z_t  # gently decreasing theta aloft
                rtm = np.full(nz, 1e-3)
                em = rng.uniform(1.0, 3.0) * np.ones(nz)
                lmax = 2000.0
            elif fam == 3:
                # saturated band -> rc_par > 0 latent-heating branch
                T = np.maximum(rng.uniform(285.0, 295.0)
                               - 0.0055 * z_t, 210.0)
                thlm = T / exner
                rsl, _ = d.drv_sat(p, T)
                rtm = 0.6 * rsl
                band = (z_t > 500.0) & (z_t < 4000.0)
                rtm[band] = 1.05 * rsl[band]
                em = rng.uniform(0.2, 0.8) * bl_m + 1e-3
            elif fam == 4:
                # bitwise-uniform dCAPE column: uniform thlm/rtm/
                # thv_ds and constant thvm => dCAPE_dz_j is IDENTICAL
                # at every level, covering the equal-dCAPE special
                # branch in both directions (sign of delta picks
                # up/down exhaustion)
                thlm = np.full(nz, 290.0)
                rtm = np.full(nz, 2e-3)
                delta = (0.05 if c % 16 < 8 else -0.05) \
                    * rng.uniform(0.5, 2.0)
                thv_ds = np.full(nz, 290.0)
                thvm = thlm + 0.61 * thv_ds * rtm + delta
                em = np.full(nz, rng.uniform(0.02, 0.3))
                for k, v in [("thvm", thvm), ("thlm", thlm),
                             ("rtm", rtm), ("em", em), ("p", p),
                             ("exner", exner), ("thv_ds", thv_ds)]:
                    ml_ins[k][c] = v
                ml_mu[c] = mu
                ml_lmax[c] = lmax
                ls, lu, ld, e = d.drv_compute_mixing_length(
                    thvm, thlm, rtm, em, lmax, mu, p, exner, thv_ds)
                assert e == 0, c
                ml_outs["lscale"][c] = ls
                ml_outs["lscale_up"][c] = lu
                ml_outs["lscale_down"][c] = ld
                continue
            elif fam == 5:
                # near-zero TKE everywhere: first-level exhaustion
                T = 288.0 + 0.006 * z_t
                thlm = T / exner
                rtm = np.full(nz, 5e-4)
                em = np.full(nz, rng.uniform(1e-4, 6e-4))
            elif fam == 6:
                # entrainment-rate / Lscale_max sweep on a convective
                # column
                T = np.maximum(300.0 - 0.0085 * z_t, 210.0)
                thlm = T / exner
                rtm = np.full(nz, 4e-3)
                em = 0.8 * bl_m + 1e-3
                mu = [0.00025, 0.001, 0.002, MU_EAMV3][c // 8 % 4]
                lmax = 1.0e5 if c % 16 >= 8 else 0.25 * HOST_DXY
            else:
                # randomized with sharp inversions + noise (nonlocal
                # Lscale_up/down smoothing coverage)
                T = np.maximum(rng.uniform(280.0, 300.0)
                               - rng.uniform(0.004, 0.009) * z_t, 205.0)
                T += 4.0 * (z_t > rng.uniform(1000.0, 4000.0))
                thlm = T / exner + rng.normal(0.0, 0.3, nz)
                rtm = np.clip(rng.uniform(0.2, 0.9)
                              * d.drv_sat(p, T)[0], 1e-7, 0.02)
                em = rng.uniform(0.05, 1.0, nz) * bl_m + 6e-4
            thv_ds = thlm.copy()
            thvm = thlm + 0.61 * thv_ds * rtm
        for k, v in [("thvm", thvm), ("thlm", thlm), ("rtm", rtm),
                     ("em", em), ("p", p), ("exner", exner),
                     ("thv_ds", thv_ds)]:
            ml_ins[k][c] = v
        ml_mu[c] = mu
        ml_lmax[c] = lmax
        ls, lu, ld, e = d.drv_compute_mixing_length(
            thvm, thlm, rtm, em, lmax, mu, p, exner, thv_ds)
        assert e == 0, c
        ml_outs["lscale"][c] = ls
        ml_outs["lscale_up"][c] = lu
        ml_outs["lscale_down"][c] = ld

    # branch coverage: Lscale_max cap and the lminh surface floor
    assert (ml_outs["lscale"] == ml_lmax[:, None]).any()
    lminh = np.maximum(0.0, 500.0 - (z_t - z_m[0])) * lmin \
        * (1.0 / 500.0)
    assert (ml_outs["lscale_up"][:, 1:] == lminh[None, 1:]).any()
    assert (ml_outs["lscale_down"][:, 1:] == lminh[None, 1:]).any()
    print(f"ml coverage: cap {(ml_outs['lscale'] == ml_lmax[:, None]).sum()}"
          f", up floor {(ml_outs['lscale_up'][:, 1:] == lminh[None, 1:]).sum()}"
          f", down floor {(ml_outs['lscale_down'][:, 1:] == lminh[None, 1:]).sum()}")

    out.update({f"ml_in_{k}": v for k, v in ml_ins.items()})
    out.update(ml_mu=ml_mu, ml_lscale_max=ml_lmax)
    out.update({f"ml_out_{k}": v for k, v in ml_outs.items()})

    # ---- BV: calc_brunt_vaisala_freq_sqd, all three formula variants
    # (dry-T0 = EAMv3, dry-thvm, moist Durran-Klemp) -------------------
    nbv = 16
    bv_ins = {k: np.zeros((nbv, nz)) for k in BV_IN}
    bv_outs = {k: np.zeros((nbv, nz)) for k in
               ["dry_t0", "dry_thvm", "moist"]}
    for c in range(nbv):
        j = c % 25
        thlm = seg_ins["thlm"][j]
        exner = seg_ins["exner"][j]
        rtm = seg_ins["rtm"][j]
        p = seg_ins["p"][j]
        thvm = seg_outs[j, :, 1]
        rcm = seg_ins["rcm"][j].copy()
        if c % 2 == 1:
            # synthetic cloud band (moist formula sensitivity)
            band = (z_t > 800.0) & (z_t < 3000.0)
            rcm[band] = np.maximum(rcm[band], 3e-4)
        for k, v in [("thlm", thlm), ("exner", exner), ("rtm", rtm),
                     ("rcm", rcm), ("p", p), ("thvm", thvm)]:
            bv_ins[k][c] = v
        args = [thlm, exner, rtm, rcm, p, thvm]
        d.drv_set_bv_flags(False, False)
        bv_outs["dry_t0"][c] = d.drv_brunt_vaisala(*args)
        d.drv_set_bv_flags(False, True)
        bv_outs["dry_thvm"][c] = d.drv_brunt_vaisala(*args)
        d.drv_set_bv_flags(True, False)
        bv_outs["moist"][c] = d.drv_brunt_vaisala(*args)
        d.drv_set_bv_flags(False, False)   # restore EAMv3
    out.update({f"bv_in_{k}": v for k, v in bv_ins.items()})
    out.update({f"bv_out_{k}": v for k, v in bv_outs.items()})

    # ---- SPLAT: term_wp2_splat / term_wp3_splat ----------------------
    # EAMv3 C_wp2_splat = 0 (exact-zero tendencies) plus a nonzero
    # sweep that exercises the five/dt clip.
    nsp = 12
    sp_ins = {k: np.zeros((nsp, nz)) for k in
              ["wp2", "wp2_zt", "wp3", "tau_zm", "tau_zt"]}
    sp_c = np.zeros(nsp)
    sp_outs = {k: np.zeros((nsp, nz)) for k in
               ["wp2_splat", "wp3_splat"]}
    for c in range(nsp):
        j = c % 25
        wp2 = seg_ins["wp2"][j].copy()
        wp3 = seg_ins["wp3"][j].copy()
        if c % 3 == 2:
            # sharp step in wp2 -> huge d(sqrt(wp2))/dz -> clip arm
            wp2 = np.where(z_m < 400.0, 1.5, W_TOL_SQD)
            wp3 = 0.5 * wp2 ** 1.5
        wp2_zt = np.maximum(_zm2zt_f(wp2, nz), W_TOL_SQD)
        tau_zm = seg_outs[j, :, 7]
        tau_zt = seg_outs[j, :, 6]
        cs = 0.0 if c < 6 else 2.0
        for k, v in [("wp2", wp2), ("wp2_zt", wp2_zt), ("wp3", wp3),
                     ("tau_zm", tau_zm), ("tau_zt", tau_zt)]:
            sp_ins[k][c] = v
        sp_c[c] = cs
        w2s, w3s = d.drv_term_splat(DT_EAMV3, cs, wp2, wp2_zt, wp3,
                                    tau_zm, tau_zt)
        sp_outs["wp2_splat"][c] = w2s
        sp_outs["wp3_splat"][c] = w3s
    # coverage: the five/dt clip must engage somewhere for C=2
    clip_lim = 5.0 / DT_EAMV3
    engaged = (sp_outs["wp2_splat"][6:] ==
               -sp_ins["wp2"][6:] * clip_lim).any()
    assert engaged
    assert (sp_outs["wp2_splat"][:6] == -0.0).all()  # C=0 exact zeros
    out.update({f"sp_in_{k}": v for k, v in sp_ins.items()})
    out.update(sp_c_wp2_splat=sp_c)
    out.update({f"sp_out_{k}": v for k, v in sp_outs.items()})

    # ---- SV: calc_surface_varnce scalar sweeps -----------------------
    sv_cases = []
    for f in SEG_SFC_FLUXES:
        sv_cases.append(dict(wpthlp=f[0], wprtp=f[1], upwp=f[2],
                             vpwp=f[3], um=5.0, vm=-3.0, lup=150.0,
                             splat=0.0, tau=600.0))
    # sign sweeps + tolerance floors + wstar branch
    for wpthlp in [-0.3, -0.02, 0.0, 1e-13, 0.02, 0.4]:
        for wprtp in [-2e-4, 0.0, 8e-4]:
            sv_cases.append(dict(wpthlp=wpthlp, wprtp=wprtp,
                                 upwp=-0.05, vpwp=0.02, um=8.0,
                                 vm=1.0, lup=300.0, splat=0.0,
                                 tau=600.0))
    # splatting correction: both arms (mild -> additive, strong ->
    # min_wp2 correlation-guard branch)
    for splat, tau in [(-1e-5, 100.0), (-0.005, 900.0), (-5.0, 900.0),
                       (0.0, 900.0)]:
        sv_cases.append(dict(wpthlp=0.15, wprtp=3e-4, upwp=-0.3,
                             vpwp=0.2, um=10.0, vm=-4.0, lup=800.0,
                             splat=splat, tau=tau))
    # randomized
    for _ in range(14):
        sv_cases.append(dict(
            wpthlp=rng.uniform(-0.4, 0.4), wprtp=rng.uniform(-1e-3, 1e-3),
            upwp=rng.uniform(-0.6, 0.6), vpwp=rng.uniform(-0.6, 0.6),
            um=rng.uniform(-15.0, 15.0), vm=rng.uniform(-15.0, 15.0),
            lup=rng.uniform(0.1, 2000.0),
            splat=-(10.0 ** rng.uniform(-6.0, 0.5)),
            tau=rng.uniform(50.0, 3000.0)))
    nsv = len(sv_cases)
    sv_in = {k: np.zeros(nsv) for k in
             ["wpthlp", "wprtp", "upwp", "vpwp", "um", "vm", "lup",
              "splat", "tau"]}
    sv_out = np.zeros((nsv, 6))
    for c, case in enumerate(sv_cases):
        for k in sv_in:
            sv_in[k][c] = case[k]
        outs, e = d.drv_calc_surface_varnce(
            case["upwp"], case["vpwp"], case["wpthlp"], case["wprtp"],
            case["um"], case["vm"], case["lup"], case["splat"],
            case["tau"])
        assert e == 0, c
        sv_out[c] = outs
    # coverage: tolerance floors, and the min_wp2 correlation guard
    assert (sv_out[:, 3] == THL_TOL ** 2).any()
    assert (sv_out[:, 4] == RT_TOL ** 2).any()
    mmcf_sqd = 0.99 * 0.99
    minv = np.maximum.reduce([
        np.full(nsv, W_TOL_SQD),
        sv_in["wprtp"] * sv_in["wprtp"] / (sv_out[:, 4] * mmcf_sqd),
        sv_in["wpthlp"] * sv_in["wpthlp"] / (sv_out[:, 3] * mmcf_sqd)])
    guard = sv_out[:, 0] == minv
    assert guard.any() and (~guard).any()
    print(f"sv coverage: thl floor {(sv_out[:, 3] == THL_TOL**2).sum()}, "
          f"rt floor {(sv_out[:, 4] == RT_TOL**2).sum()}, "
          f"min_wp2 guard {guard.sum()}/{nsv}")
    out.update({f"sv_in_{k}": v for k, v in sv_in.items()})
    out.update(sv_outs=sv_out)

    return out


def main(which=None):
    params, idx = eamv3_params()
    xp2_idx = {n: int(i) - 1 for n, i in
               zip(XP2_IDX_NAMES, d.drv_param_indices_xp2())}
    lscale_idx = {n: int(i) - 1 for n, i in
                  zip(LSCALE_IDX_NAMES, d.drv_param_indices_lscale())}
    sha = subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"],
                         capture_output=True, text=True).stdout.strip()
    meta = dict(
        source_sha=sha,
        fflags="-O2 -fPIC -ffree-line-length-none "
               "-fallow-argument-mismatch -std=legacy -ffp-contract=off "
               "-DCLUBB_CAM -DCLUBB_SGS -DCLUBB_REAL_TYPE=dp",
        lapack="Ubuntu reference liblapack (container)",
        eamv3_overrides=EAMV3_OVERRIDES,
        param_indices={k: int(v) for k, v in idx.items()},
        xp2_param_indices=xp2_idx,
        lscale_param_indices=lscale_idx,
        lscale_config=dict(
            # model_flags defaults; clubb_intr never overrides these
            # (EAMv3 clubb_stabcorrect=F only touches
            # l_stability_correct_Kh_N2_zm / l_diffuse_rtm_and_thlm)
            l_stability_correct_tau_zm=True,
            l_diag_Lscale_from_tau=False,
            l_use_C7_Richardson=False, l_use_C11_Richardson=False,
            l_use_wp3_pr3=False,          # => Cx_fnc_Richardson = 0
            l_avg_Lscale=False,           # COMPILE-TIME parameter in
                                          # advance_clubb_core: the
                                          # perturbed-Lscale calls and
                                          # averaging are dead in EAM
            l_Lscale_plume_centered=False, l_use_ice_latent=False,
            l_brunt_vaisala_freq_moist=False,
            l_use_thvm_in_bv_freq=False, l_sat_mixrat_lookup=False,
            l_include_ice=False,          # compute_rsat_parcel
            l_andre_1978=False,           # calc_surface_varnce
            l_tke_aniso=True, lmin=4.0, T0=300.0,
            lscale_max=0.25 * HOST_DXY, mu=MU_EAMV3),
        seg_slots=SEG_SLOTS, sfc_slots=SFC_SLOTS,
        xp2_config=dict(l_iter_xp2_xpyp=True, l_single_C2_Skw=False,
                        l_explicit_turbulent_adv_xpyp=False,
                        l_upwind_xpyp_ta=True,
                        l_min_xp2_from_corr_wx=False,
                        l_C2_cloud_frac=False, l_hole_fill=True,
                        l_tke_aniso=True, l_up2_vp2_sponge_damp=False,
                        l_clip_large_rtp2=True, rtp2_clip_coef=0.5,
                        gamma_over_implicit_ts=1.5),
        config=dict(grid_type=3, l_implemented=True, sclr_dim=0,
                    hydromet_dim=0, edsclr_dim=29, theta0=300.0,
                    ts_nudge=86400.0, saturation_formula="flatau",
                    iiPDF_type="ADG1", l_stats=False, debug_level=0,
                    dt=DT_EAMV3, ipdf_call_placement=2,
                    l_do_expldiff_rtm_thlm=True,
                    l_trapezoidal_rule_zt=True, l_trapezoidal_rule_zm=True,
                    l_call_pdf_closure_twice=True, l_use_cloud_cover=True,
                    l_use_ice_latent=False, l_rcm_supersat_adj=False,
                    l_rtm_nudge=False, l_gamma_Skw=True),
        outs_slots=OUTS_SLOTS,
    )
    meta_s = json.dumps(meta)

    gens = {
        "grid": lambda: gen_grid(),
        "tridag": lambda: gen_tridag(),
        "sat": lambda: gen_sat(),
        "pdf_closure": lambda: gen_pdf_closure(params, idx),
        "pdf_driver": lambda: gen_pdf_driver(params, idx),
        "xp2": lambda: gen_xp2(params, idx),
        "lscale": lambda: gen_lscale(params, idx),
    }
    for name in (which or gens):
        np.savez_compressed(GOLDEN / f"clubb_{name}.npz", meta=meta_s,
                            params=params, **gens[name]())
        print(f"wrote clubb_{name}.npz")


if __name__ == "__main__":
    import sys
    main(sys.argv[1:] or None)
