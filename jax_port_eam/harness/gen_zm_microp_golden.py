#!/usr/bin/env python3
"""Tier-1 golden generator for the ZM convective microphysics
(zm_microphysics.F90 zm_mphy + zm_conv_main with zm_microp=.true.;
run in the scream-dev container).

Two golden sets in one archive:

config k -- direct zm_mphy kernel on 40 synthetic gathered plume
  columns sweeping every regime family (5 families x 8 columns):
    0-7   warm plume        (whole updraft above 0C: warm rain only)
    8-15  mixed-phase plume (freezing/riming/graupel active)
    16-23 ice plume         (cold surface, cmei-dominated)
    24-31 aerosol sweep     (na_scale 0.05 .. 20 on a mixed plume)
    32-39 dynamics sweep    (weak..strong updrafts, shallow/deep,
                             cond 5e-8..2e-6, one cond=0 no-op column,
                             one eps0=0 no-op column, one cmei-gap
                             column for the libase boundary branch)

config m/n -- full zm_conv_main with zm_param%zm_microp=.true. on the
  same 42-column synthetic sounding set as gen_zm_conv_golden.py
  (identical recipe + rng seed), EAMv3 phys="default" params plus the
  microphysics tuning (zmconv_auto_fac=7.0, zmconv_accr_fac=1.5,
  zmconv_micro_dcs=150e-6; old_snow=.false. as zm_conv_intr sets it
  under zmconv_microp), with column-varying MAM4 aerosol profiles.
  m: is_first_step=T (CAPE trigger); n: is_first_step=F (DCAPE).

Aerosol activation strategy (PORTING_PLAN.md row 9): the REAL
activate_drop_mam + nucleate_ice_conv modules are compiled into the
harness; only bulk ndrop_bam is an abort stub (never executed: the
modal scheme is EAMv3 production). The MAM4 mode/species data below
(sigmag, densities, hygroscopicities, dgnum lo/hi) stand in for
rad_constituents physprop data and are recorded in the metadata --
they are pure inputs on both sides of the comparison.

deltat (zm_mphy activation timescale factors) = 1800 s via the
time_manager stub. All indices stored 1-based as Fortran returns
them. state/diag/microp field orders are documented in
drivers/zm_microp_driver.F90 and repeated in the metadata.
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "fmod"))
import eam_zm_microp_f  # noqa: E402

d = eam_zm_microp_f.zm_microp_driver

GRAV, CP, RD = 9.80616, 1004.64, 287.042
TIME_STEP = 1800.0
nlev = 72
F = lambda a: np.asarray(a, dtype=np.float64, order="F")  # noqa: E731

# ---- MAM4 (no-MOM CPP variant) aerosol configuration ----
NMODES, NSPECMX = 4, 6
SIGMAG = np.array([1.8, 1.6, 1.8, 1.6])          # accum, aitken, coarse, pcarbon
NSPEC = np.array([6, 3, 3, 2], dtype=np.int32)
M_ACCUM, M_AITKEN, M_COARSE = 1, 2, 3            # 1-based mode indices
L_DUST, L_NACL, L_SO4 = 1, 2, 3                  # coarse-mode species indices
SPECDENS = np.zeros((NSPECMX, NMODES))
SPECHYGRO = np.zeros((NSPECMX, NMODES))
SPECDENS[:, 0] = [1770., 1000., 1000., 1700., 2600., 1900.]   # so4,pom,soa,bc,dst,ncl
SPECHYGRO[:, 0] = [0.507, 1e-10, 0.14, 1e-10, 0.068, 1.16]
SPECDENS[:3, 1] = [1770., 1000., 1900.]                       # so4,soa,ncl
SPECHYGRO[:3, 1] = [0.507, 0.14, 1.16]
SPECDENS[:3, 2] = [2600., 1900., 1770.]                       # dst,ncl,so4
SPECHYGRO[:3, 2] = [0.068, 1.16, 0.507]
SPECDENS[:2, 3] = [1000., 1700.]                              # pom,bc
SPECHYGRO[:2, 3] = [1e-10, 1e-10]
DGNUMLO = np.array([5.35e-8, 8.7e-9, 1.0e-6, 1.0e-8])
DGNUMHI = np.array([4.4e-7, 5.2e-8, 4.0e-6, 1.0e-7])
# zm_aero_init: voltonumblo from dgnumLO (upper number bound),
# voltonumbhi from dgnumHI (lower number bound)
V2NLO = 1.0 / (np.pi / 6.0 * DGNUMLO**3 * np.exp(4.5 * np.log(SIGMAG)**2))
V2NHI = 1.0 / (np.pi / 6.0 * DGNUMHI**3 * np.exp(4.5 * np.log(SIGMAG)**2))
DGNUM = np.array([1.1e-7, 2.6e-8, 2.0e-6, 5.0e-8])
NUM0 = np.array([6.0e8, 4.0e8, 2.0e5, 1.0e8])    # [#/kg] baseline
MMR0 = np.zeros((NSPECMX, NMODES))
MMR0[:6, 0] = [8e-10, 4e-10, 6e-10, 5e-11, 2e-10, 3e-10]
MMR0[:3, 1] = [6e-11, 3e-11, 1e-11]
MMR0[:3, 2] = [2e-9, 1.2e-9, 8e-11]
MMR0[:2, 3] = [2e-10, 4e-11]

AUTO_FAC, ACCR_FAC, MICRO_DCS = 7.0, 1.5, 150.0e-6
MUCON, DCON = 5.3, 25.0e-6
LAMBDADPCU0, MUDPCU0 = (MUCON + 1.0) / DCON, MUCON

d.drv_init(int(TIME_STEP), SIGMAG)

STATE_FIELDS = ["qc", "qi", "nc", "ni", "qcde", "qide", "qnide",
                "ncde", "nide", "nsde", "qni", "qr", "ns", "nr", "qg",
                "ng", "rprd", "sprd", "frz", "wu", "lamc", "pgam"]
DIAG_FIELDS = ["autolm", "accrlm", "bergnm", "fhtimm", "fhtctm",
               "fhmlm", "hmpim", "accslm", "dlfm", "autoln", "accrln",
               "bergnn", "fhtimn", "fhtctn", "fhmln", "accsln",
               "activn", "dlfn", "autoim", "accsim", "difm", "nuclin",
               "autoin", "accsin", "hmpin", "difn", "trspcm", "trspcn",
               "trspim", "trspin", "accgrm", "accglm", "accgslm",
               "accgsrm", "accgirm", "accgrim", "accgrsm", "accgsln",
               "accgsrn", "accgirn", "accsrim", "acciglm", "accigrm",
               "accsirm", "accigln", "accigrn", "accsirn", "accgln",
               "accgrn", "accilm", "acciln", "fallrm", "fallsm",
               "fallgm", "fallrn", "fallsn", "fallgn", "fhmrm", "dsfm",
               "dsfn"]
MICROP_FIELDS = ["wu", "qliq", "qice", "qrain", "qsnow", "qgraupel",
                 "qnl", "qni", "qnr", "qns", "qng", "sprd", "mudpcu",
                 "lambdadpcu", "qcde", "qide", "qsde", "ncde", "nide",
                 "nsde", "dif", "dsf", "dnlf", "dnif", "dnsf", "frz",
                 "cmel", "cmei", "autolm", "accrlm", "bergnm",
                 "fhtimm", "fhtctm", "fhmlm", "hmpim", "accslm",
                 "dlfm", "dsfm", "autoln", "accrln", "bergnn",
                 "fhtimn", "fhtctn", "fhmln", "accsln", "activn",
                 "dlfn", "dsfn", "autoim", "accsim", "difm", "nuclin",
                 "autoin", "accsin", "hmpin", "difn", "trspcm",
                 "trspcn", "trspim", "trspin", "accgrm", "accglm",
                 "accgslm", "accgsrm", "accgirm", "accgrim", "accgrsm",
                 "accgsln", "accgsrn", "accgirn", "accsrim", "acciglm",
                 "accigrm", "accsirm", "accigln", "accigrn", "accsirn",
                 "accgln", "accgrn", "accilm", "acciln", "fallrm",
                 "fallsm", "fallgm", "fallrn", "fallsn", "fallgn",
                 "fhmrm"]

out = {}

# ===================================================================
# config k: direct zm_mphy kernel golden
# ===================================================================
rngk = np.random.default_rng(20260713)
nk = 40
MSG_K = 2

t_srf = np.empty(nk)
cond = np.empty(nk)
na_scale = np.ones(nk)
su_off = np.empty(nk)
jb = np.empty(nk, dtype=np.int32)
jt = np.empty(nk, dtype=np.int32)
jlcl = np.empty(nk, dtype=np.int32)
eps0 = np.full(nk, 1.0e-4)
cmei_gap = np.zeros(nk, dtype=bool)

fam = np.repeat(np.arange(5), 8)
for i in range(nk):
    f, j = fam[i], i % 8
    if f == 0:    # warm plume: shallow, updraft never below 0C
        t_srf[i], su_off[i] = 305.0 + j * 0.5, 1.5
        cond[i] = 2.0e-7 + 4.0e-8 * j
        jb[i], jt[i], jlcl[i] = 63, 48 - j, 59
    elif f == 1:  # mixed-phase deep plume
        t_srf[i], su_off[i] = 300.0 + 0.5 * j, 1.0 + 0.3 * j
        cond[i] = 3.0e-7 + 6.0e-8 * j
        jb[i], jt[i], jlcl[i] = 61, 22 + j, 55 - (j % 3)
    elif f == 2:  # ice plume
        t_srf[i], su_off[i] = 268.0 + 0.8 * j, 1.0 + 0.2 * j
        cond[i] = 1.5e-7 + 4.0e-8 * j
        jb[i], jt[i], jlcl[i] = 62, 26 + j, 57
    elif f == 3:  # aerosol sweep on a mixed plume
        t_srf[i], su_off[i] = 301.0, 1.8
        cond[i] = 4.0e-7
        na_scale[i] = [0.05, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0][j]
        jb[i], jt[i], jlcl[i] = 61, 25, 55
    else:         # dynamics sweep
        t_srf[i] = 301.0
        su_off[i] = [0.2, 0.5, 1.0, 2.0, 4.0, 2.0, 2.0, 2.0][j]
        cond[i] = [4e-7, 4e-7, 4e-7, 4e-7, 4e-7, 5e-8, 2e-6, 0.0][j]
        jb[i], jt[i], jlcl[i] = 61, 25, 55
        if j == 5:
            jb[i], jt[i], jlcl[i] = 58, 44, 55   # shallow
        if j == 4:
            eps0[i] = 0.0                        # no-op column
        if j == 6:
            cmei_gap[i] = True                   # libase branch

zi1 = np.linspace(1.0, 0.0, nlev + 1) ** 1.4 * 16000.0
zf = np.broadcast_to(zi1, (nk, nlev + 1)).copy()
zf *= 1 + 0.02 * rngk.uniform(-1, 1, (nk, 1))
zm = 0.5 * (zf[:, :-1] + zf[:, 1:])
te = np.maximum(t_srf[:, None] - 6.5e-3 * zm, 200.0) \
    + rngk.uniform(-0.3, 0.3, (nk, nlev))
pm = 1000.0 * np.exp(-zm / 7600.0)
qe = 0.016 * np.exp(-zm / 2200.0) * (1 + 0.05 * rngk.uniform(-1, 1, (nk, nlev)))
su = te + GRAV / CP * zf[:, :nlev] + su_off[:, None]
qu = qe * 1.05
kk1 = np.arange(1, nlev + 1)[None, :]
inpl = (kk1 > jt[:, None]) & (kk1 <= jb[:, None])
mu_top = 0.15 + 0.5 * rngk.uniform(0, 1, nk)
mu = np.where(kk1 >= jb[:, None], 1.0,
              np.where(inpl, mu_top[:, None]
                       + (1.0 - mu_top[:, None]) * (kk1 - jt[:, None])
                       / (jb - jt)[:, None], 0.0))
dz = zf[:, :nlev] - zf[:, 1:]
dmu = np.zeros((nk, nlev))
dmu[:, :-1] = mu[:, :-1] - mu[:, 1:]
eu = np.where(inpl & (dmu > 0), dmu / dz, 0.0)
du = np.where(inpl & (dmu < 0), -dmu / dz, 0.0)
rows = np.arange(nk)
du[rows, jt - 1] = mu[rows, jt] / dz[rows, jt - 1]
gamhat = np.clip(0.622 * 2.5e6**2 * qe / (RD * te**2 * CP), 0.05, 4.0)
tu = su - GRAV / CP * zf[:, :nlev]
fice_k = np.clip((273.15 - tu) / 40.0, 0.0, 1.0)
cu_k = np.where(inpl, cond[:, None], 0.0)
cmel = cu_k * (1 - fice_k)
cmei = cu_k * fice_k
# cmei-gap column: interrupt the ice condensation mid-plume so the
# libase (cmei(k-1)>qsmall & cmei(k)<qsmall) boundary branch fires
gap = cmei_gap[:, None] & (kk1 >= 38) & (kk1 <= 42)
cmei = np.where(gap, 0.0, cmei)

zprof = np.exp(-zm / 2500.0)
numg = np.zeros((nk, nlev, NMODES))
mmrg = np.zeros((nk, nlev, NSPECMX, NMODES))
for m in range(NMODES):
    prof = zprof if m != M_COARSE - 1 else np.exp(-zm / 1800.0)
    numg[:, :, m] = NUM0[m] * na_scale[:, None] * (0.15 + prof)
    for l in range(NSPEC[m]):
        mmrg[:, :, l, m] = MMR0[l, m] * na_scale[:, None] * (0.15 + prof)
dgn = np.broadcast_to(DGNUM, (nk, nlev, NMODES)).copy()
dgn *= 1 + 0.1 * rngk.uniform(-1, 1, (nk, nlev, NMODES))

state, diag = d.drv_zm_mphy(
    MSG_K, AUTO_FAC, ACCR_FAC, MICRO_DCS, LAMBDADPCU0, MUDPCU0,
    jb, jt, jlcl, F(su), F(qu), F(mu), F(du), F(eu), F(zf), F(pm),
    F(te), F(qe), F(gamhat), eps0, F(cmel), F(cmei),
    NSPEC, M_ACCUM, M_AITKEN, M_COARSE, L_DUST, L_NACL, L_SO4,
    SIGMAG[1], F(SPECDENS), F(SPECHYGRO), V2NLO, V2NHI,
    F(numg), F(mmrg), F(dgn))
state = np.asarray(state)
diag = np.asarray(diag)

sk = dict(zip(STATE_FIELDS, np.moveaxis(state, 2, 0)))
# sanity: positivity, regime limits, aerosol monotonicity
for nm in ["qc", "qi", "qr", "qni", "qg", "nc", "ni", "nr", "ns", "ng"]:
    assert np.all(np.isfinite(sk[nm])) and np.all(sk[nm] >= 0.0), nm
warmf = slice(0, 8)
assert sk["qi"][warmf].max() == 0.0 and sk["sprd"][warmf].max() == 0.0, \
    "warm family must be liquid-only"
assert sk["rprd"][warmf].max() > 0.0
mixedf = slice(8, 16)
assert sk["qi"][mixedf].max() > 0.0 and sk["sprd"][mixedf].max() > 0.0
icef = slice(16, 24)
# fully-glaciated-from-base plumes never fire the kqi ice boundary
# condition (cmei>qsmall at cloud base already), so qi stays 0 while
# the freezing heating (frz = cmei + ...) is still active -- a
# faithful regime worth goldening
assert sk["frz"][icef].max() > 0.0
aer = slice(24, 32)
ncmax = sk["nc"][aer].max(axis=1)
assert ncmax[-1] > ncmax[0], "more aerosol -> more droplets"
i_eps0 = 32 + 4
i_cond0 = 32 + 7
assert np.all(state[i_eps0, :, :20] == 0.0), "eps0=0 must be a no-op"
assert np.all(state[i_cond0, :, :20] == 0.0), "cond=0 must be a no-op"

out.update(k_state=state, k_diag=diag, k_jb=jb, k_jt=jt, k_jlcl=jlcl,
           k_su=su, k_qu=qu, k_mu=mu, k_du=du, k_eu=eu, k_zf=zf,
           k_pm=pm, k_te=te, k_qe=qe, k_gamhat=gamhat, k_eps0=eps0,
           k_cmel=cmel, k_cmei=cmei, k_numg=numg, k_mmrg=mmrg,
           k_dgnum=dgn, k_msg=np.array(MSG_K), k_family=fam)
print(f"config k: max qc={sk['qc'].max():.3e} qi={sk['qi'].max():.3e} "
      f"rprd={sk['rprd'].max():.3e} sprd={sk['sprd'].max():.3e} "
      f"wu={sk['wu'].max():.2f}")

# ===================================================================
# configs m/n: full zm_conv_main with zm_microp=.true.
# (sounding recipe identical to gen_zm_conv_golden.py, seed included)
# ===================================================================
rng = np.random.default_rng(20260712)
nfam, npc = 6, 7
ncol = nfam * npc

ai = np.linspace(0.0, 1.0, nlev + 1) ** 1.7
pint1 = 225.5 + ai * (1.0e5 - 225.5)
pint = np.broadcast_to(pint1, (ncol, nlev + 1)).copy()
pint *= (1 + 0.03 * rng.uniform(-1, 1, (ncol, 1)))
pmid = 0.5 * (pint[:, :-1] + pint[:, 1:])
pdel = pint[:, 1:] - pint[:, :-1]
ps = pmid[:, -1:]
limcnv = int(np.argmax(pint1 >= 4000.0))


def qs_bolton(T, p):
    es = 611.2 * np.exp(17.67 * (T - 273.15) / (T - 29.65))
    return 0.622 * es / (p - es)


T = np.zeros((ncol, nlev))
rh = np.zeros((ncol, nlev))
sig = pmid / ps
for f, (ts, ex, rh0) in enumerate([
        (303.0, 0.19, 0.85), (296.0, 0.16, 0.70), (282.0, 0.10, 0.50),
        (291.0, 0.18, 0.80), (312.0, 0.21, 0.05), (300.0, 0.19, 0.99)]):
    rowsl = slice(f * npc, (f + 1) * npc)
    tsv = ts + rng.uniform(-1.0, 2.0, (npc, 1))
    T[rowsl] = np.maximum(tsv * sig[rowsl] ** ex, 195.0)
    rh[rowsl] = np.clip(rh0 * sig[rowsl] ** 0.8, 0.02, rh0)
T += rng.uniform(-0.5, 0.5, T.shape)
q = np.maximum(rh * qs_bolton(T, pmid), 1.0e-9)
q *= 1 + 0.05 * rng.uniform(-1, 1, q.shape)

tv = T * (1 + 0.608 * q)
zint = np.zeros((ncol, nlev + 1))
zmid = np.zeros((ncol, nlev))
for k in range(nlev - 1, -1, -1):
    zmid[:, k] = zint[:, k + 1] + RD * tv[:, k] / GRAV \
        * np.log(pint[:, k + 1] / pmid[:, k])
    zint[:, k] = zint[:, k + 1] + RD * tv[:, k] / GRAV \
        * np.log(pint[:, k + 1] / pint[:, k])
z_srf = rng.uniform(0.0, 500.0, (ncol,))
geos = GRAV * z_srf
pbl_hgt = np.where(np.arange(ncol) % 3 == 0, 2200.0, 1000.0)
tpert = np.where(np.arange(ncol) % 2 == 0, 0.0,
                 rng.uniform(0.2, 1.5, ncol))
landfrac = np.tile([0.0, 1.0, 0.37, 1.0, 0.0, 1.0, 0.5],
                   nfam)[:ncol].astype(np.float64)
omega = rng.uniform(-0.6, 0.6, (ncol, nlev))
t_star = T + rng.uniform(-0.5, 0.5, T.shape)
q_star = np.maximum(q * (1 + 0.03 * rng.uniform(-1, 1, q.shape)), 1.0e-9)

# column-varying MAM4 aerosol: clean maritime .. polluted continental,
# with a dust-heavy subset for contact freezing
na_col = np.tile([0.1, 0.3, 1.0, 2.0, 5.0, 10.0, 1.0], nfam)[:ncol]
dust_col = np.tile([1.0, 1.0, 0.2, 5.0, 1.0, 0.5, 20.0], nfam)[:ncol]
zprof_m = np.exp(-zmid / 2500.0)
num_a = np.zeros((ncol, nlev, NMODES))
mmr_a = np.zeros((ncol, nlev, NSPECMX, NMODES))
for m in range(NMODES):
    scale = na_col if m != M_COARSE - 1 else na_col * dust_col
    prof = zprof_m if m != M_COARSE - 1 else np.exp(-zmid / 1800.0)
    num_a[:, :, m] = NUM0[m] * scale[:, None] * (0.15 + prof)
    for l in range(NSPEC[m]):
        mmr_a[:, :, l, m] = MMR0[l, m] * scale[:, None] * (0.15 + prof)
dgn_m = np.broadcast_to(DGNUM, (ncol, nlev, NMODES)).copy()
dgn_m *= 1 + 0.1 * rng.uniform(-1, 1, dgn_m.shape)

P_V3 = dict(tau=3600.0, alfa=0.14, ke=2.5e-6, dmpdz=-0.7e-3,
            tpert_fix=1, tpert_fac=2.0, tiedke_add=0.8,
            c0_lnd=0.0020, c0_ocn=0.0020, num_cin=1,
            mx_bot_lyr_adj=1, trig_dcape=1, trig_ull=1,
            clos_dyn_adj=1, old_snow=0, no_deep_pbl=0,
            auto_fac=AUTO_FAC, accr_fac=ACCR_FAC, micro_dcs=MICRO_DCS)
NAMES = ("lengath", "gather_index", "msemax_klev_g", "jctop", "jcbot",
         "jt", "prec", "heat", "qtnd", "cape", "dcape", "mcon", "pflx",
         "zdu", "mflx_up", "entr_up", "detr_up", "mflx_dn", "entr_dn",
         "p_del", "dsubcld", "ql", "rliq", "rprd", "dlf",
         "microp", "rice")

out.update(t=T, q=q, omega=omega, pmid=pmid, pint=pint, pdel=pdel,
           geos=geos, zmid=zmid, zint=zint, pbl_hgt=pbl_hgt,
           tpert=tpert, landfrac=landfrac, t_star=t_star,
           q_star=q_star, limcnv=np.array(limcnv),
           time_step=np.array(TIME_STEP),
           num_a=num_a, mmr_a=mmr_a, dgnum=dgn_m,
           aero_sigmag=SIGMAG, aero_nspec=NSPEC,
           aero_specdens=SPECDENS, aero_spechygro=SPECHYGRO,
           aero_voltonumblo=V2NLO, aero_voltonumbhi=V2NHI)

for cfg, first in (("m", 1), ("n", 0)):
    p = P_V3
    res = d.drv_zm_conv_main_microp(
        TIME_STEP, first, limcnv, p["no_deep_pbl"],
        p["tau"], p["alfa"], p["ke"], p["dmpdz"], p["tpert_fix"],
        p["tpert_fac"], p["tiedke_add"], p["c0_lnd"], p["c0_ocn"],
        p["num_cin"], p["mx_bot_lyr_adj"], p["trig_dcape"],
        p["trig_ull"], p["clos_dyn_adj"], p["old_snow"],
        p["auto_fac"], p["accr_fac"], p["micro_dcs"],
        F(T), F(q), F(omega), F(pmid), F(pint), F(pdel), geos,
        F(zmid), F(zint), pbl_hgt, tpert, landfrac,
        F(t_star), F(q_star),
        NSPEC, M_ACCUM, M_AITKEN, M_COARSE, L_DUST, L_NACL, L_SO4,
        SIGMAG[1], F(SPECDENS), F(SPECHYGRO), V2NLO, V2NHI,
        F(num_a), F(mmr_a), F(dgn_m))
    r = dict(zip(NAMES, res))
    out.update({f"{n}_{cfg}": np.asarray(v) for n, v in r.items()})
    ng = int(r["lengath"])
    mp = np.asarray(r["microp"])
    sprd = mp[:, :, MICROP_FIELDS.index("sprd")]
    dif = mp[:, :, MICROP_FIELDS.index("dif")]
    print(f"cfg {cfg}: lengath={ng} "
          f"gathered={np.sort(r['gather_index'][:ng])}")
    print(f"  prec[mm/day]={np.array2string(r['prec'] * 8.64e7, precision=1)}")
    print(f"  sprd_max={sprd.max():.3e} dif_max={dif.max():.3e} "
          f"rice_max={r['rice'].max():.3e}")
    assert np.all(np.isfinite(mp))
    assert ng > 0 and r["prec"].max() > 0.0
    # snow/ice production somewhere in the deep tropical columns
    assert sprd.max() > 0.0 and dif.max() > 0.0

sha = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                     text=True, cwd="/work/E3SM").stdout.strip()
meta = {"scheme": "zm_microphysics zm_mphy kernel (config k) + "
                  "zm_conv_main with zm_microp=.true. (configs m, n)",
        "source_sha": sha, "nlev": nlev, "ncol_kernel": nk,
        "ncol_main": ncol, "time_step": TIME_STEP, "index_base": 1,
        "limcnv": limcnv, "msg_main": limcnv - 1,
        "msg_kernel": MSG_K,
        "activation": "real activate_drop_mam + nucleate_ice_conv "
                      "compiled (modal scheme); ndrop_bam abort stub "
                      "(bulk never executed)",
        "cpp_variant": "no MODAL_AERO_4MODE_MOM / RAIN_EVAP_TO_COARSE_"
                       "AERO: wght = dmc/(ssmc+dmc+so4mc)",
        "erf_gamma": "gfortran: shr_spfn erf/gamma resolve to glibc "
                     "intrinsics (HAVE_*_INTRINSICS via __GFORTRAN__)",
        "params": {**P_V3, "limcnv": limcnv,
                   "lambdadpcu0": LAMBDADPCU0, "mudpcu0": MUDPCU0,
                   "deltat": TIME_STEP},
        "aero": {"nmodes": NMODES, "nspecmx": NSPECMX,
                 "sigmag": SIGMAG.tolist(), "nspec": NSPEC.tolist(),
                 "mode_idx": [M_ACCUM, M_AITKEN, M_COARSE],
                 "coarse_species_idx": [L_DUST, L_NACL, L_SO4],
                 "dgnumlo": DGNUMLO.tolist(), "dgnumhi": DGNUMHI.tolist(),
                 "dgnum": DGNUM.tolist(),
                 "note": "MAM4 physprop stand-ins; pure data inputs "
                         "recorded here and replayed identically"},
        "state_fields": STATE_FIELDS, "diag_fields": DIAG_FIELDS,
        "microp_fields": MICROP_FIELDS,
        "families_kernel": ["warm", "mixed", "ice", "aerosol_sweep",
                            "dynamics_sweep"]}
gold = Path(__file__).resolve().parents[1] / "golden" / "zm_microp_golden.npz"
np.savez_compressed(gold, __metadata__=json.dumps(meta), **out)
print(f"wrote {gold}")
