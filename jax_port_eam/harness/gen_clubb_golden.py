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

Run `python3 gen_clubb_golden.py [archive ...]` with archive names
(grid/tridag/sat/pdf_closure/pdf_driver) to regenerate a subset;
no arguments regenerates everything.

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
# (bld/namelist_files/namelist_defaults_eam.xml; clubb_wpxp_L_thresh
# is the general default at line 871).  clubb_C2rt=1.75 triggers the
# read_parameters coupling C2thl=C2rt, C2rtthl=1.3*C2rt.
EAMV3_OVERRIDES = {
    "C1": 2.4, "C1b": 2.8, "C1c": 0.75,
    "C2rt": 1.75, "C2thl": 1.75, "C2rtthl": 1.75 * 1.3,
    "C6rtb": 7.50, "C6rtc": 0.50, "C6thlb": 7.50, "C6thlc": 0.50,
    "C8": 5.2, "C11": 0.70, "C11b": 0.20, "C11c": 0.85,
    "gamma_coef": 0.12, "gamma_coefb": 0.28, "gamma_coefc": 1.2,
    "mu": 0.0005, "wpxp_L_thresh": 60.0,
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


def main(which=None):
    params, idx = eamv3_params()
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
    }
    for name in (which or gens):
        np.savez_compressed(GOLDEN / f"clubb_{name}.npz", meta=meta_s,
                            params=params, **gens[name]())
        print(f"wrote clubb_{name}.npz")


if __name__ == "__main__":
    import sys
    main(sys.argv[1:] or None)
