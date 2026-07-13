#!/usr/bin/env python3
"""Container smoke test for eam_zm_microp_f: wrapper dtype check +
hand-checkable physics on the direct zm_mphy kernel (warm-rain
autoconversion grows with condensation supply, freezing only at/below
0C, mass/number positivity) and one full zm_conv_main microp=T call."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "fmod"))
import eam_zm_microp_f  # noqa: E402

d = eam_zm_microp_f.zm_microp_driver

# f2py kind-map check: r8 must be double everywhere
assert "array('d')" in d.drv_zm_mphy.__doc__, d.drv_zm_mphy.__doc__
assert "array('d')" in d.drv_zm_conv_main_microp.__doc__

# ---- MAM4 (no-MOM variant) aerosol configuration ----
NMODES, NSPECMX = 4, 6
SIGMAG = np.array([1.8, 1.6, 1.8, 1.6])
NSPEC = np.array([6, 3, 3, 2], dtype=np.int32)
M_ACCUM, M_AITKEN, M_COARSE = 1, 2, 3
L_DUST, L_NACL, L_SO4 = 1, 2, 3
SPECDENS = np.zeros((NSPECMX, NMODES))
SPECHYGRO = np.zeros((NSPECMX, NMODES))
# accum: so4, pom, soa, bc, dst, ncl
SPECDENS[:, 0] = [1770., 1000., 1000., 1700., 2600., 1900.]
SPECHYGRO[:, 0] = [0.507, 1e-10, 0.14, 1e-10, 0.068, 1.16]
# aitken: so4, soa, ncl
SPECDENS[:3, 1] = [1770., 1000., 1900.]
SPECHYGRO[:3, 1] = [0.507, 0.14, 1.16]
# coarse: dst, ncl, so4
SPECDENS[:3, 2] = [2600., 1900., 1770.]
SPECHYGRO[:3, 2] = [0.068, 1.16, 0.507]
# pcarbon: pom, bc
SPECDENS[:2, 3] = [1000., 1700.]
SPECHYGRO[:2, 3] = [1e-10, 1e-10]
DGNUMLO = np.array([5.35e-8, 8.7e-9, 1.0e-6, 1.0e-8])
DGNUMHI = np.array([4.4e-7, 5.2e-8, 4.0e-6, 1.0e-7])
# zm_aero_init: voltonumblo from dgnumLO (upper number bound)
V2NLO = 1.0 / (np.pi / 6.0 * DGNUMLO**3 * np.exp(4.5 * np.log(SIGMAG)**2))
V2NHI = 1.0 / (np.pi / 6.0 * DGNUMHI**3 * np.exp(4.5 * np.log(SIGMAG)**2))
DGNUM = np.array([1.1e-7, 2.6e-8, 2.0e-6, 5.0e-8])

d.drv_init(1800, SIGMAG)

GRAV, CP, RD = 9.80616, 1004.64, 287.042
nlev = 72
F = lambda a: np.asarray(a, dtype=np.float64, order="F")  # noqa: E731


def make_plume(t_srf=302.0, mu_top=0.2, cond=4.0e-7, jb=61, jt=26,
               jlcl=56, na_scale=1.0):
    """One synthetic gathered column of zm_mphy inputs (1-based idx)."""
    ncol = 1
    zi = np.linspace(1.0, 0.0, nlev + 1) ** 1.4 * 16000.0  # interfaces, top->sfc
    zf = zi[None, :].copy()
    # hydrostatic-ish environment
    te = np.maximum(t_srf - 6.5e-3 * 0.5 * (zi[:-1] + zi[1:]), 200.0)[None, :]
    pm = (1000.0 * np.exp(-0.5 * (zi[:-1] + zi[1:]) / 7600.0))[None, :]
    qe = 0.016 * np.exp(-0.5 * (zi[:-1] + zi[1:]) / 2200.0)[None, :]
    su = te + GRAV / CP * zf[:, :nlev] + 1.0   # updraft ~1K warmer
    qu = qe * 1.05
    k = np.arange(1, nlev + 1)
    inplume = (k > jt) & (k <= jb)
    mu = np.where(k >= jb, 1.0,
                  np.where(inplume, mu_top + (1.0 - mu_top)
                           * (k - jt) / (jb - jt), 0.0))[None, :]
    mu[0, :jt - 1] = 0.0
    dz = zf[0, :nlev] - zf[0, 1:]
    du = np.zeros((1, nlev))
    eu = np.zeros((1, nlev))
    # crude flux-consistent entrain/detrain
    dmu = np.zeros(nlev)
    dmu[:-1] = mu[0, :-1] - mu[0, 1:]
    eu[0] = np.where(inplume & (dmu > 0), dmu / dz, 0.0)
    du[0] = np.where(inplume & (dmu < 0), -dmu / dz, 0.0)
    du[0, jt - 1] = mu[0, jt] / dz[jt - 1]
    gamhat = np.clip(0.622 * 2.5e6**2 * qe / (RD * te**2 * CP), 0.05, 4.0)
    eps0 = np.array([1.0e-4])
    # condensation: liquid where updraft T>0C, ice below
    tu = su - GRAV / CP * zf[:, :nlev]
    fice = np.clip((273.15 - tu) / 40.0, 0.0, 1.0)
    cu = np.where(inplume, cond, 0.0)[None, :]
    cmel = cu * (1 - fice)
    cmei = cu * fice
    numg = np.zeros((ncol, nlev, NMODES))
    numg[:, :, 0] = 6.0e8 * na_scale
    numg[:, :, 1] = 4.0e8 * na_scale
    numg[:, :, 2] = 2.0e5 * na_scale
    numg[:, :, 3] = 1.0e8 * na_scale
    mmrg = np.zeros((ncol, nlev, NSPECMX, NMODES))
    mmrg[:, :, :NSPEC[0], 0] = np.array([8e-10, 4e-10, 6e-10, 5e-11,
                                         2e-10, 3e-10]) * na_scale
    mmrg[:, :, :NSPEC[1], 1] = np.array([6e-11, 3e-11, 1e-11]) * na_scale
    mmrg[:, :, :NSPEC[2], 2] = np.array([2e-9, 1.2e-9, 8e-11]) * na_scale
    mmrg[:, :, :NSPEC[3], 3] = np.array([2e-10, 4e-11]) * na_scale
    dgn = np.broadcast_to(DGNUM, (ncol, nlev, NMODES)).copy()
    args = (2, 7.0, 1.5, 150.0e-6, (5.3 + 1.0) / 25.0e-6, 5.3,
            np.array([jb], dtype=np.int32), np.array([jt], dtype=np.int32),
            np.array([jlcl], dtype=np.int32),
            F(su), F(qu), F(mu), F(du), F(eu), F(zf), F(pm), F(te),
            F(qe), F(gamhat), eps0, F(cmel), F(cmei),
            NSPEC, M_ACCUM, M_AITKEN, M_COARSE, L_DUST, L_NACL, L_SO4,
            SIGMAG[1], F(SPECDENS), F(SPECHYGRO), V2NLO, V2NHI,
            F(numg), F(mmrg), F(dgn))
    return args


SNAMES = ["qc", "qi", "nc", "ni", "qcde", "qide", "qnide", "ncde",
          "nide", "nsde", "qni", "qr", "ns", "nr", "qg", "ng", "rprd",
          "sprd", "frz", "wu", "lamc", "pgam"]

st1, dg1 = d.drv_zm_mphy(*make_plume(cond=2.0e-7))
st2, dg2 = d.drv_zm_mphy(*make_plume(cond=8.0e-7))
s1 = dict(zip(SNAMES, np.moveaxis(st1, 2, 0)))
s2 = dict(zip(SNAMES, np.moveaxis(st2, 2, 0)))

for nm in ["qc", "qi", "nc", "ni", "qr", "qni", "qg", "nr", "ns", "ng"]:
    assert np.all(np.isfinite(s1[nm])) and np.all(s1[nm] >= 0.0), nm
assert s1["qc"].max() > 0.0, "expected in-plume cloud water"
assert s2["rprd"].sum() > s1["rprd"].sum(), \
    "precip production should grow with condensation supply"
# freezing only where the updraft is at/below 0C: frz is stored at
# k-1 for plume level k, so allow one level of offset
tu = make_plume()[9] - GRAV / CP * make_plume()[14][:, :nlev]
warm = tu[0] > 274.5
assert np.all(np.abs(s2["frz"][0, warm][1:]) < 1e-30), \
    "no freezing in warm layers"
cold_plume = s2["qi"][0] > 0
assert s2["frz"].max() > 0.0, "expected freezing in the cold plume"
print("zm_mphy smoke OK:",
      f"max qc={s2['qc'].max():.3e} max qi={s2['qi'].max():.3e}",
      f"max rprd={s2['rprd'].max():.3e} max sprd={s2['sprd'].max():.3e}",
      f"max wu={s2['wu'].max():.2f}")

# ---- full zm_conv_main with zm_microp=.true. on one warm column ----
rng = np.random.default_rng(7)
ncol = 4
ai = np.linspace(0.0, 1.0, nlev + 1) ** 1.7
pint = (225.5 + ai * (1.0e5 - 225.5))[None, :].repeat(ncol, 0)
pmid = 0.5 * (pint[:, :-1] + pint[:, 1:])
pdel = pint[:, 1:] - pint[:, :-1]
sig = pmid / pmid[:, -1:]
T = np.maximum(303.0 * sig ** 0.19, 195.0) + rng.uniform(-0.3, 0.3, pmid.shape)
es = 611.2 * np.exp(17.67 * (T - 273.15) / (T - 29.65))
q = np.maximum(0.85 * np.clip(0.85 * sig ** 0.8, 0.02, 0.85)
               * 0.622 * es / (pmid - es), 1e-9)
tv = T * (1 + 0.608 * q)
zint = np.zeros((ncol, nlev + 1))
zmid = np.zeros((ncol, nlev))
for k in range(nlev - 1, -1, -1):
    zmid[:, k] = zint[:, k + 1] + RD * tv[:, k] / GRAV \
        * np.log(pint[:, k + 1] / pmid[:, k])
    zint[:, k] = zint[:, k + 1] + RD * tv[:, k] / GRAV \
        * np.log(pint[:, k + 1] / pint[:, k])
limcnv = int(np.argmax(pint[0] >= 4000.0))
geos = np.zeros(ncol)
num_a = np.zeros((ncol, nlev, NMODES))
num_a[:, :, :] = np.array([6.0e8, 4.0e8, 2.0e5, 1.0e8])
mmr_a = np.zeros((ncol, nlev, NSPECMX, NMODES))
mmr_a[:, :, :NSPEC[0], 0] = [8e-10, 4e-10, 6e-10, 5e-11, 2e-10, 3e-10]
mmr_a[:, :, :NSPEC[1], 1] = [6e-11, 3e-11, 1e-11]
mmr_a[:, :, :NSPEC[2], 2] = [2e-9, 1.2e-9, 8e-11]
mmr_a[:, :, :NSPEC[3], 3] = [2e-10, 4e-11]
dgn = np.broadcast_to(DGNUM, (ncol, nlev, NMODES)).copy()

res = d.drv_zm_conv_main_microp(
    1800.0, 1, limcnv, 0, 3600.0, 0.14, 2.5e-6, -0.7e-3, 1, 2.0, 0.8,
    0.0020, 0.0020, 1, 1, 1, 1, 1, 0, 7.0, 1.5, 150.0e-6,
    F(T), F(q), F(np.zeros_like(T)), F(pmid), F(pint), F(pdel), geos,
    F(zmid), F(zint), np.full(ncol, 1000.0), np.zeros(ncol),
    np.zeros(ncol), F(T), F(q),
    NSPEC, M_ACCUM, M_AITKEN, M_COARSE, L_DUST, L_NACL, L_SO4,
    SIGMAG[1], F(SPECDENS), F(SPECHYGRO), V2NLO, V2NHI,
    F(num_a), F(mmr_a), F(dgn))
lengath = int(res[0])
prec, cape = res[6], res[9]
microp, rice = res[-2], res[-1]
assert lengath > 0, "expected trigger"
assert np.all(np.isfinite(microp)) and np.all(np.isfinite(prec))
assert prec.max() > 0.0
sprd = microp[:, :, 11]
print(f"zm_conv_main microp OK: lengath={lengath} "
      f"cape_max={cape.max():.0f} prec_max_mmday={prec.max() * 8.64e7:.2f} "
      f"sprd_max={sprd.max():.3e} rice={rice}")
