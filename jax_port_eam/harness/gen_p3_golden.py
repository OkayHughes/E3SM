#!/usr/bin/env python3
"""Tier-1 golden generator for the EAM P3 stratiform microphysics
(micro_p3.F90 p3_main, eam variant; run in the scream-dev container).

40 columns x 72 levels, five regime families of 8:
  0- 7 warm rain            (liquid cloud 1-4 km + drizzle, T>0C where
                             hydrometeors live; autoconv/accretion/evap/
                             cloud+rain sedimentation)
  8-15 mixed-phase          (cloud spanning 0..-25C: riming, hetero
                             freezing (CNT or Bigg), Bergeron, melting)
 16-23 ice-only             (cold columns, ice supersaturated:
                             deposition, nucleation, ice sed, aggregation)
 24-31 heavy precip         (deep cloud, qr up to ~4e-3, ice aloft,
                             melting through 0C, strong sedimentation)
 32-39 evaporation/edge     (sub-saturated rain shafts; sub-qsmall
                             wisps -> clipping; one fully clear column;
                             one cold clear ice-supersat column
                             (nucleation-only); one T>+2C ice column
                             (ice_complete_melting); one T<-40C liquid
                             column (homogeneous freezing))

Configs (all chained: each step feeds P3's output state + the
interface's t_prev/qv_prev = end-of-step T/qv back in):
  a: use_hetfrz_classnuc=T, do_predict_nc=T   3 steps  (EAMv3 default)
  b: use_hetfrz_classnuc=F, do_predict_nc=T   2 steps  (Bigg immersion)
  c: use_hetfrz_classnuc=T, do_predict_nc=F   1 step   (Cooper nuc + nccnst)
  d: hetfrz=T, predict_nc=T, do_Cooper_inP3=T 1 step   (Cooper add-on)
  e: hetfrz=T, do_prescribed_CCN=T            1 step   (prescribed CCN)

Parameters: EAMv3 phys="default" (namelist_defaults_eam.xml):
p3_autocon_coeff=30500, p3_qc_autocon_expon=3.19,
p3_nc_autocon_expon=-1.10, p3_accret_coeff=117.25,
p3_qc_accret_expon=1.15, p3_wbf_coeff=1.0, p3_mincdnc=20e6,
p3_max_mean_rain_size=0.005, p3_embryonic_rain_size=25e-6,
micro_nccons=200e6, do_precip_off=F.

Table: p3_lookup_table_1.dat-v4.1.1 (the only version on local disk;
EAMv3 default 4.1.2 lives in atm/cam/physprops, not present locally).
The table is pure input data recorded in this archive, so the version
cancels between Fortran and JAX. Rain tables from p3_get_tables
(p3_init_b with mu_r_constant = 0 — EAM differs from SCREAM's 1.0
here) are stored for exact replay.
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "fmod"))
import eam_p3_f  # noqa: E402

d = eam_p3_f.p3_driver

TABLE_DIR = "/work/e3sm-inputdata/atm/scream/tables"
TABLE_VERSION = "4.1.1"

# physconst values (physconst stub, verbatim shr_const derivations)
RAIR = 6.02214e26 * 1.38065e-23 / 28.966
CP, GRAV = 1.00464e3, 9.80616
QSMALL, MINCLD = 1.0e-14, 0.0001
DT = 1800.0
ncol, nlev = 40, 72
F = lambda a: np.asarray(a, dtype=np.float64, order="F")  # noqa: E731

P_V3 = dict(autocon_coeff=30500.0, accret_coeff=117.25,
            qc_autocon_expon=3.19, nc_autocon_expon=-1.10,
            qc_accret_expon=1.15, wbf_coeff=1.0, mincdnc=20.0e6,
            max_mean_rain_size=0.005, embryonic_rain_size=25.0e-6,
            nccnst=200.0e6)

STATE_FIELDS = ["qc", "nc", "qr", "nr", "th", "qv", "qi", "qm", "ni", "bm"]
DIAG_FIELDS = ["diag_eff_radius_qc", "diag_eff_radius_qi", "rho_qi",
               "qv2qi_depos_tend", "precip_total_tend", "nevapr",
               "qr_evap_tend", "mu_c", "lamc", "liq_ice_exchange",
               "vap_liq_exchange", "vap_ice_exchange",
               "diag_equiv_reflectivity", "diag_ze_rain", "diag_ze_ice"]
FLUX_FIELDS = ["precip_liq_flux", "precip_ice_flux", "rflx", "sflx", "cflx"]

d.drv_init(TABLE_DIR.ljust(256), TABLE_VERSION.ljust(16), 1, 0)
mu_r_tab, revap_tab, vn_tab, vm_tab = (np.asarray(a)
                                       for a in d.drv_get_tables())

# ------------------------------------------------------------------
# qv_sat probe grid (unit bisection of the saturation port)
# ------------------------------------------------------------------
t_probe = np.linspace(170.0, 330.0, 161)
p_probe = np.geomspace(2.0e3, 1.05e5, 161)
qsat_liq = np.asarray(d.drv_qv_sat(t_probe, p_probe, 0))
qsat_ice = np.asarray(d.drv_qv_sat(t_probe, p_probe, 1))

# ------------------------------------------------------------------
# synthetic soundings: 5 families x 8 columns
# ------------------------------------------------------------------
rng = np.random.default_rng(20260712)

zi1 = np.linspace(1.0, 0.0, nlev + 1) ** 1.35 * 24000.0
zi = np.broadcast_to(zi1, (ncol, nlev + 1)).copy()
zi *= 1 + 0.02 * rng.uniform(-1, 1, (ncol, 1))
zm = 0.5 * (zi[:, :-1] + zi[:, 1:])
dz = zi[:, :-1] - zi[:, 1:]

fam = np.repeat(np.arange(5), 8)
t_srf = np.empty(ncol)
rh0 = np.empty(ncol)
for i in range(ncol):
    f, j = fam[i], i % 8
    if f == 0:
        t_srf[i], rh0[i] = 296.0 + j, 0.85
    elif f == 1:
        t_srf[i], rh0[i] = 283.0 + 0.8 * j, 0.80
    elif f == 2:
        t_srf[i], rh0[i] = 258.0 + 1.2 * j, 0.75
    elif f == 3:
        t_srf[i], rh0[i] = 299.0 + 0.5 * j, 0.90
    else:
        t_srf[i], rh0[i] = 292.0 + j, 0.45

lapse = 6.5e-3 + 0.3e-3 * rng.uniform(-1, 1, (ncol, 1))
T = np.maximum(t_srf[:, None] - lapse * zm, 200.0) \
    + rng.uniform(-0.3, 0.3, (ncol, nlev))
pm = 1.0e5 * np.exp(-zm / 7600.0)


def qsat_bolton(T, p):
    es = 611.2 * np.exp(17.67 * (T - 273.15) / (T - 29.65))
    return 0.622 * es / np.maximum(p - es, 1.0)


rh_prof = np.clip(rh0[:, None] * (1.0 - 0.55 * zm / 24000.0), 0.03, 0.99)
qv = np.minimum(rh_prof * qsat_bolton(T, pm), 0.030)

exner = 1.0 / ((pm * 1e-5) ** (RAIR / CP))     # EAM "exner": th = T*exner
th = T * exner
rho = pm / (RAIR * T)
dpres = rho * GRAV * dz                        # p3 part1: rho = dpres/dz/g

# hydrometeor placement per family
qc = np.zeros((ncol, nlev)); nc = np.zeros((ncol, nlev))
qr = np.zeros((ncol, nlev)); nr = np.zeros((ncol, nlev))
qi = np.zeros((ncol, nlev)); ni = np.zeros((ncol, nlev))
qm = np.zeros((ncol, nlev)); bm = np.zeros((ncol, nlev))

for i in range(ncol):
    f, j = fam[i], i % 8
    Ti = T[i]
    if f == 0:      # warm cloud 1-4 km, drizzle below 2.5 km
        lay = (zm[i] > 1000) & (zm[i] < 4000) & (Ti > 273.15)
        qc[i][lay] = (2e-4 + 1e-4 * j) * (1 + 0.1 * rng.uniform(-1, 1, lay.sum()))
        nc[i][lay] = 3e7 + 1e7 * j
        shaft = (zm[i] < 2500) & (Ti > 273.15)
        qr[i][shaft] = 2e-5 + 2e-5 * j
        nr[i][shaft] = 2e4 + 1e4 * j
    elif f == 1:    # mixed-phase 2-7 km
        lay = (zm[i] > 2000) & (zm[i] < 7000)
        qc[i][lay] = 1.5e-4 + 0.6e-4 * j
        nc[i][lay] = 5e7
        ice = lay & (Ti < 271.0)
        qi[i][ice] = 0.8e-4 + 0.5e-4 * j
        ni[i][ice] = 3e4 + 2e4 * j
        qm[i][ice] = (0.2 + 0.07 * j) * qi[i][ice]
        bm[i][ice] = qm[i][ice] / (300.0 + 30.0 * j)
        shaft = zm[i] < 3500
        qr[i][shaft] = 3e-5 + 1e-5 * j
        nr[i][shaft] = 3e4
    elif f == 2:    # ice-only 3-10 km
        lay = (zm[i] > 3000) & (zm[i] < 10000)
        qi[i][lay] = 1e-4 + 0.8e-4 * j
        ni[i][lay] = 2e4 + 3e4 * j
        qm[i][lay] = (0.1 + 0.1 * j) * qi[i][lay]
        bm[i][lay] = qm[i][lay] / (250.0 + 50.0 * j)
        qv[i][lay] = np.minimum(1.12 * qsat_bolton(Ti[lay], pm[i][lay])
                                * 0.95, 0.03)  # ice supersaturated
    elif f == 3:    # heavy precip: deep cloud + big qr + ice aloft
        lay = (zm[i] > 800) & (zm[i] < 9000)
        qc[i][lay] = 6e-4 + 2e-4 * j
        nc[i][lay] = 1e8
        shaft = zm[i] < 6000
        qr[i][shaft] = 8e-4 + 4e-4 * j
        nr[i][shaft] = 3e5 + 2e5 * j
        ice = lay & (Ti < 268.0)
        qi[i][ice] = 4e-4 + 2e-4 * j
        ni[i][ice] = 8e4
        qm[i][ice] = 0.45 * qi[i][ice]
        bm[i][ice] = qm[i][ice] / 420.0
    else:           # evap / edge cases
        if j == 0:
            pass  # fully clear
        elif j == 1:
            # cold clear ice-supersaturated: nucleation-only path
            lay = (zm[i] > 6000) & (zm[i] < 11000)
            qv[i][lay] = np.minimum(1.10 * qsat_bolton(Ti[lay], pm[i][lay])
                                    * 0.95, 0.03)
        elif j == 2:
            # sub-qsmall wisps: part1 clipping branches
            lay = (zm[i] > 2000) & (zm[i] < 5000)
            qc[i][lay] = 3e-15
            nc[i][lay] = 1e2
            qi[i][lay & (Ti < 270)] = 5e-15
            qr[i][lay] = 2e-15
        elif j == 3:
            # tiny warm ice (qsmall < qi < 1e-8, T > 0C): melt-to-rain clip
            lay = (zm[i] > 500) & (zm[i] < 2000) & (Ti > 274.0)
            qi[i][lay] = 5e-9
            ni[i][lay] = 1e3
        elif j == 4:
            # ice at T > +2C: ice_complete_melting
            lay = (zm[i] > 300) & (zm[i] < 1500) & (Ti > 276.0)
            qi[i][lay] = 2e-4
            ni[i][lay] = 4e4
            qm[i][lay] = 0.05 * qi[i][lay]  # lightly rimed
            bm[i][lay] = qm[i][lay] / 300.0
        elif j == 5:
            # supercooled liquid below -40C: homogeneous freezing
            lay = (Ti < 230.0) & (zm[i] < 12000)
            qc[i][lay] = 1e-4
            nc[i][lay] = 4e7
            qr[i][lay] = 2e-5
            nr[i][lay] = 2e4
        else:
            # sub-saturated rain shafts (strong evaporation)
            shaft = (zm[i] < 4000)
            qr[i][shaft] = (1e-4 + 1e-4 * (j - 6))
            nr[i][shaft] = 1e5
            # dry column: rh0 = 0.45 family default already

# interface cloud fractions (get_cloud_fraction, max_overlap)
ast = np.clip((qc + qi) * 1.5e3, 0.0, 0.9)
ast = np.where((qc + qi) > 1e-9, np.maximum(ast, 0.15), ast)
cld_frac_l = np.maximum(ast, MINCLD)
cld_frac_i = np.maximum(ast, MINCLD)
cld_frac_r = np.maximum(ast, MINCLD)
for k in range(1, nlev):
    src = (qr[:, k - 1] >= QSMALL) | (qi[:, k - 1] >= QSMALL)
    cld_frac_r[:, k] = np.where(src, np.maximum(cld_frac_r[:, k - 1],
                                                cld_frac_r[:, k]),
                                cld_frac_r[:, k])

# aerosol/activation inputs
zprof = np.exp(-zm / 3000.0)
nc_nuceat_tend = 3.0e4 * zprof * (0.2 + rng.uniform(0, 1, (ncol, 1)))
ni_activated = 8.0e4 * np.exp(-((zm - 8000.0) / 3500.0) ** 2) \
    * (0.1 + rng.uniform(0, 1, (ncol, 1)))
frzimm = 0.02 * zprof * (0.1 + rng.uniform(0, 1, (ncol, 1)))
frzcnt = 0.002 * zprof * (0.1 + rng.uniform(0, 1, (ncol, 1)))
frzdep = 0.001 * zprof * (0.1 + rng.uniform(0, 1, (ncol, 1)))
inv_qc_relvar = np.clip(rng.uniform(0.5, 8.0, (ncol, nlev)), 0.001, 10.0)
nccn_prescribed = 6.0e7 * zprof * (0.3 + rng.uniform(0, 1, (ncol, 1)))
col_location = np.stack([np.arange(1.0, ncol + 1),
                         rng.uniform(-180, 180, ncol),
                         rng.uniform(-90, 90, ncol)], axis=1)

state0 = {"qc": qc, "nc": nc, "qr": qr, "nr": nr, "th": th, "qv": qv,
          "qi": qi, "qm": qm, "ni": ni, "bm": bm}

out = {
    "zm": zm, "dz": dz, "pres": pm, "dpres": dpres, "exner": exner,
    "cld_frac_l": cld_frac_l, "cld_frac_i": cld_frac_i,
    "cld_frac_r": cld_frac_r, "nc_nuceat_tend": nc_nuceat_tend,
    "nccn_prescribed": nccn_prescribed, "ni_activated": ni_activated,
    "frzimm": frzimm, "frzcnt": frzcnt, "frzdep": frzdep,
    "inv_qc_relvar": inv_qc_relvar, "col_location": col_location,
    "t_init": T, "family": fam,
    "mu_r_table": mu_r_tab, "revap_table": revap_tab,
    "vn_table": vn_tab, "vm_table": vm_tab,
    "qsat_probe_t": t_probe, "qsat_probe_p": p_probe,
    "qsat_probe_liq": qsat_liq, "qsat_probe_ice": qsat_ice,
}


def run_config(cfg, hetfrz, predict_nc, prescribed_ccn, do_cooper, nsteps):
    d.drv_set_flags(hetfrz, do_cooper)
    st = {k: v.copy() for k, v in state0.items()}
    t_prev = T.copy()
    qv_prev = st["qv"].copy()
    for s in range(nsteps):
        out[f"{cfg}_in_state_{s}"] = np.stack(
            [st[k] for k in STATE_FIELDS], axis=-1)
        out[f"{cfg}_in_tprev_{s}"] = t_prev.copy()
        out[f"{cfg}_in_qvprev_{s}"] = qv_prev.copy()
        state, diag, flux, tend, surf = d.drv_p3_main(
            DT, s + 1, predict_nc, prescribed_ccn, 0,
            P_V3["autocon_coeff"], P_V3["accret_coeff"],
            P_V3["qc_autocon_expon"], P_V3["nc_autocon_expon"],
            P_V3["qc_accret_expon"], P_V3["wbf_coeff"], P_V3["mincdnc"],
            P_V3["max_mean_rain_size"], P_V3["embryonic_rain_size"],
            P_V3["nccnst"],
            F(st["qc"]), F(st["nc"]), F(st["qr"]), F(st["nr"]),
            F(st["th"]), F(st["qv"]), F(st["qi"]), F(st["qm"]),
            F(st["ni"]), F(st["bm"]),
            F(pm), F(dz), F(nc_nuceat_tend), F(nccn_prescribed),
            F(ni_activated), F(frzimm), F(frzcnt), F(frzdep),
            F(inv_qc_relvar), F(dpres), F(exner),
            F(cld_frac_r), F(cld_frac_l), F(cld_frac_i),
            F(qv_prev), F(t_prev), F(col_location))
        state = np.asarray(state); diag = np.asarray(diag)
        flux = np.asarray(flux); tend = np.asarray(tend)
        surf = np.asarray(surf)
        assert np.all(np.isfinite(state)), f"{cfg} step {s}: non-finite state"
        out[f"{cfg}_state_{s}"] = state
        out[f"{cfg}_diag_{s}"] = diag
        out[f"{cfg}_flux_{s}"] = flux
        out[f"{cfg}_tend_{s}"] = tend
        out[f"{cfg}_surf_{s}"] = surf
        for k, j in zip(STATE_FIELDS, range(10)):
            st[k] = state[:, :, j].copy()
        # micro_p3_tend: t_prev/qv_prev = end-of-step T (th/exner) and qv
        t_prev = st["th"] / exner
        qv_prev = st["qv"].copy()
    sk = dict(zip(STATE_FIELDS, np.moveaxis(out[f"{cfg}_state_0"], 2, 0)))
    prl = out[f"{cfg}_surf_{nsteps-1}"][:, 0]
    pri = out[f"{cfg}_surf_{nsteps-1}"][:, 1]
    print(f"cfg {cfg}: prec_liq[mm/d] max={prl.max()*8.64e7:.2f} "
          f"prec_ice max={pri.max()*8.64e7:.2f} "
          f"qi_max={sk['qi'].max():.2e} qr_max={sk['qr'].max():.2e}")
    return prl, pri


prl_a, pri_a = run_config("a", 1, 1, 0, 0, 3)
run_config("b", 0, 1, 0, 0, 2)
run_config("c", 1, 0, 0, 0, 1)
run_config("d", 1, 1, 0, 1, 1)
run_config("e", 1, 1, 1, 0, 1)

# ------------------------------------------------------------------
# sedimentation sub-kernel goldens (drv_cloud/rain/ice_sed).
# The full-chain state comparison is limited by knife-edge CFL
# substep amplification of sub-1e-11 process noise, so the sed
# kernels get their own tight golden: inputs derived from the initial
# state (pure data; both sides receive them identically).
# ------------------------------------------------------------------
rho_sed = dpres / dz / GRAV
inv_rho_sed = 1.0 / rho_sed
RHO_1000MB = 1.0e5 / (RAIR * 273.15)
RHO_600MB = 6.0e4 / (RAIR * 253.15)
rhofacr = (RHO_1000MB * inv_rho_sed) ** 0.54
rhofaci = (RHO_600MB * inv_rho_sed) ** 0.54
mu_air = 1.496e-6 * T ** 1.5 / (T + 120.0)
acn = GRAV * 1000.0 / (18.0 * mu_air)
inv_dz = 1.0 / dz
zeros = np.zeros((ncol, nlev))

sed_in = dict(sed_rho=rho_sed, sed_inv_rho=inv_rho_sed,
              sed_rhofacr=rhofacr, sed_rhofaci=rhofaci, sed_acn=acn,
              sed_qc_incld=qc / cld_frac_l, sed_nc_incld=nc / cld_frac_l,
              sed_qr_incld=qr / cld_frac_r, sed_nr_incld=nr / cld_frac_r,
              sed_qi_incld=qi / cld_frac_i, sed_ni_incld=ni / cld_frac_i,
              sed_qm_incld=qm / cld_frac_i, sed_bm_incld=bm / cld_frac_i)
out.update(sed_in)

res = d.drv_cloud_sed(DT, 1, F(sed_in["sed_qc_incld"]), F(rho_sed),
                      F(inv_rho_sed), F(cld_frac_l), F(acn), F(inv_dz),
                      F(qc), F(nc), F(sed_in["sed_nc_incld"]),
                      F(zeros), F(zeros))
for name, v in zip(["qc", "nc", "nc_incld", "mu_c", "lamc",
                    "precip_liq_surf", "cflx", "qc_tend", "nc_tend"], res):
    out[f"csed_{name}"] = np.asarray(v)

res = d.drv_rain_sed(DT, P_V3["max_mean_rain_size"],
                     F(sed_in["sed_qr_incld"]), F(rho_sed),
                     F(inv_rho_sed), F(rhofacr), F(cld_frac_r), F(inv_dz),
                     F(qr), F(nr), F(sed_in["sed_nr_incld"]),
                     F(zeros), F(zeros))
for name, v in zip(["qr", "nr", "nr_incld", "mu_r", "lamr",
                    "precip_liq_surf", "precip_liq_flux", "rflx",
                    "qr_tend", "nr_tend"], res):
    out[f"rsed_{name}"] = np.asarray(v)

res = d.drv_ice_sed(DT, F(rho_sed), F(inv_rho_sed), F(rhofaci),
                    F(cld_frac_i), F(inv_dz), F(qi),
                    F(sed_in["sed_qi_incld"]), F(ni),
                    F(sed_in["sed_ni_incld"]), F(qm),
                    F(sed_in["sed_qm_incld"]), F(bm),
                    F(sed_in["sed_bm_incld"]))
for name, v in zip(["qi", "ni", "qm", "bm", "qi_incld", "ni_incld",
                    "qm_incld", "bm_incld", "precip_ice_surf", "sflx",
                    "qi_tend", "ni_tend"], res):
    out[f"ised_{name}"] = np.asarray(v)
print("sed sub-kernel goldens written")

# regime sanity on config a
warm = slice(0, 8)
sa = out["a_state_0"]
# levels above the complete-melting threshold (+2C) stay ice-free; cold
# aloft may nucleate cirrus (cell-average process needing no seed ice)
# and ice may sediment one level past 0C before melting catches it
warm_lay = T[warm] > 275.15
assert sa[warm, :, STATE_FIELDS.index("qi")][warm_lay].max() == 0.0, \
    "warm family T>+2C levels must stay ice-free"
assert pri_a[warm].max() == 0.0, "warm family must not produce ice precip"
assert prl_a[warm].max() > 0.0, "warm family must rain"
assert prl_a[24:32].max() > prl_a[:8].max(), "heavy family must rain hardest"
assert pri_a[8:24].max() > 0.0, "mixed/ice families must produce ice precip"
i_clear = 32
assert prl_a[i_clear] == 0.0 and pri_a[i_clear] == 0.0
assert np.array_equal(sa[i_clear, :, STATE_FIELDS.index("qc")],
                      state0["qc"][i_clear]), "clear col qc must be untouched"
# nucleation-only column must gain ice
i_nuc = 33
assert out["a_state_0"][i_nuc, :, STATE_FIELDS.index("qi")].max() > 0.0

sha = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                     text=True, cwd="/work/E3SM").stdout.strip()
meta = {
    "scheme": "micro_p3.F90 p3_main (eam variant)",
    "source_sha": sha, "ncol": ncol, "nlev": nlev, "dt": DT,
    "table_dir": TABLE_DIR, "table_version": TABLE_VERSION,
    "table_note": "EAMv3 default is 4.1.2 (atm/cam/physprops), not on "
                  "local disk; 4.1.1 used as pure input data recorded "
                  "here and read identically by both sides",
    "params": P_V3,
    "configs": {
        "a": {"use_hetfrz_classnuc": True, "do_predict_nc": True,
              "do_prescribed_ccn": False, "do_cooper": False, "nsteps": 3},
        "b": {"use_hetfrz_classnuc": False, "do_predict_nc": True,
              "do_prescribed_ccn": False, "do_cooper": False, "nsteps": 2},
        "c": {"use_hetfrz_classnuc": True, "do_predict_nc": False,
              "do_prescribed_ccn": False, "do_cooper": False, "nsteps": 1},
        "d": {"use_hetfrz_classnuc": True, "do_predict_nc": True,
              "do_prescribed_ccn": False, "do_cooper": True, "nsteps": 1},
        "e": {"use_hetfrz_classnuc": True, "do_predict_nc": True,
              "do_prescribed_ccn": True, "do_cooper": False, "nsteps": 1},
    },
    "state_fields": STATE_FIELDS, "diag_fields": DIAG_FIELDS,
    "flux_fields": FLUX_FIELDS,
    "families": ["warm", "mixed", "ice", "heavy", "evap_edge"],
    "compiler": "gfortran -O2 -ffp-contract=off (fbuild defaults)",
}
gold = Path(__file__).resolve().parents[1] / "golden" / "p3_golden.npz"
np.savez_compressed(gold, __metadata__=json.dumps(meta), **out)
print(f"wrote {gold}")
