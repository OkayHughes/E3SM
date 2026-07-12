#!/usr/bin/env python3
"""Tier-1 golden generator for the ZM convective transport routines
(zm_transport_tracer / zm_transport_momentum, run in the scream-dev
container).

Mass fluxes are PHYSICALLY CONSISTENT: the existing eam_zm_conv_f
extension is re-run on the zm_conv golden profiles (config a: EAMv3
phys="default" params, is_first_step=T) and its gathered outputs
mu/md/du/eu/ed/dp/jt/mx/ideep feed the transport exactly as
zm_conv_intr.F90 passes them, so the transport sees real plume
structures with triggered (lengath=18) and untriggered (24, incl. the
whole stable family) columns. The regenerated fluxes are asserted
bitwise-equal to the zm_conv golden archive and stored here so the
replay test needs only this file.

Tracers (ncnst=7; m=1 is water vapor, always skipped by the scheme):
  m=2 boundary-layer tracer (exp decay with height), fracis=1
  m=3 mid-tropospheric bump, partially soluble (fracis in [0.3,1])
  m=4 stratospheric tracer increasing with height, with small
      NEGATIVE patches in the lowest levels (minc<0 -> arithmetic
      chat branch), fracis=1
  m=5 uniform 1e-6 (transport no-op by plume mass continuity)
  m=6 'dry' mixing-ratio type (cnst_get_type_byind branch; gathered
      dpdry = 0.985*dp as state%pdeldry/100 would be)
  m=7 near-zero (1e-14) inside the convective layers, large below --
      drives tendencies that would go negative, to exercise the
      zm_microp conservation fixer

Configurations:
  t1: zm_microp=F, doconvtran=[F,T,T,T,T,T,T], dt=1800
  t2: zm_microp=T (negative-tracer fixer ON), same tracers/flags,
      dt=3600 (ztodt = 2*delta_t; the min(chat,const) flux limiters
      keep positive tracers above -0.65x/1800s, so the larger legal
      step is what makes the fixer actually fire on m=7)
  t3: zm_microp=F, doconvtran=[F,T,F,T,F,T,F], dt=900 (skip branch)
  mom1: sheared jet u + turning v, dt=1800
  mom2: uniform winds (7,-3) -- exact no-op on gathered columns
zm_transport_momentum has no runtime switch: momcu=momcd=0.4 are
compile-time parameters in zm_transport.F90 (no zmconv_mom*
namelist), nwind=2 always.
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "fmod"))
import eam_zm_conv_f  # noqa: E402
import eam_zm_transport_f  # noqa: E402

dtt = eam_zm_transport_f.zm_transport_driver.drv_transport_tracer
dtm = eam_zm_transport_f.zm_transport_driver.drv_transport_momentum
assert "array('d')" in dtt.__doc__ and "array('d')" in dtm.__doc__

F = lambda a: np.asarray(a, dtype=np.float64, order="F")  # noqa: E731

# --- re-run zm_conv (config a) for the plume mass fluxes ---
g = np.load(ROOT / "golden" / "zm_conv_golden.npz")
zmeta = json.loads(str(g["__metadata__"]))
p = zmeta["params"]["a"]
d = eam_zm_conv_f.zm_conv_driver
d.drv_init()
res = d.drv_zm_conv_main(
    float(g["time_step"]), p["is_first_step"], p["limcnv"],
    p["no_deep_pbl"], p["tau"], p["alfa"], p["ke"], p["dmpdz"],
    p["tpert_fix"], p["tpert_fac"], p["tiedke_add"], p["c0_lnd"],
    p["c0_ocn"], p["num_cin"], p["mx_bot_lyr_adj"], p["trig_dcape"],
    p["trig_ull"], p["clos_dyn_adj"], p["old_snow"],
    F(g["t"]), F(g["q"]), F(g["omega"]), F(g["pmid"]), F(g["pint"]),
    F(g["pdel"]), g["geos"], F(g["zmid"]), F(g["zint"]),
    g["pbl_hgt"], g["tpert"], g["landfrac"], F(g["t_star"]),
    F(g["q_star"]))
(lengath, ideep, mx_g, _jctop, _jcbot, jt_g, _prec, _heat, _qtnd,
 _cape, _dcape, _mcon, _pflx, _zdu, mu, eu, du, md, ed, dp, _dsub,
 _ql, _rliq, _rprd, _dlf) = res
lengath = int(lengath)
ncol, nlev = np.asarray(g["t"]).shape
print(f"zm_conv rerun: lengath={lengath} of ncol={ncol}")
assert 0 < lengath < ncol
# tie to the zm_conv golden archive (same .so + inputs => bitwise)
for name, arr in (("mflx_up", mu), ("entr_up", eu), ("detr_up", du),
                  ("mflx_dn", md), ("entr_dn", ed), ("p_del", dp),
                  ("jt", jt_g), ("msemax_klev_g", mx_g),
                  ("gather_index", ideep)):
    assert np.array_equal(np.asarray(arr), g[f"{name}_a"]), name

# --- synthetic tracers ---
rng = np.random.default_rng(20260713)
ncnst = 7
sig = np.asarray(g["pmid"]) / np.asarray(g["pmid"])[:, -1:]
q = np.zeros((ncol, nlev, ncnst))
q[:, :, 0] = np.asarray(g["q"])                       # water vapor
q[:, :, 1] = 1.0e-6 * np.exp(4.0 * (sig - 1.0)) \
    * (1.0 + 0.3 * rng.uniform(-1, 1, (ncol, nlev)))  # BL tracer
q[:, :, 2] = 1.0e-9 + 2.0e-6 * np.exp(-0.5 * ((sig - 0.45) / 0.12)**2)
q[:, :, 3] = 5.0e-7 * np.exp(-6.0 * sig) + 1.0e-10
q[:, -6:, 3] -= 2.5e-9 * rng.uniform(0.5, 1.5, (ncol, 6))  # negatives
q[:, :, 4] = 1.0e-6                                   # uniform
q[:, :, 5] = 3.0e-7 * (1.0 + np.cos(6.0 * sig)) \
    * (1.0 + 0.2 * rng.uniform(-1, 1, (ncol, nlev))) + 1.0e-12
q[:, :, 6] = np.where(sig < 0.85, 1.0e-14, 5.0e-7)    # fixer bait
assert (q[:, -6:, 3] < 0).any()

fracis = np.ones((ncol, nlev, ncnst))
fracis[:, :, 2] = 0.3 + 0.7 * sig
dpdry = np.zeros((ncol, nlev))
dpdry[:lengath] = 0.985 * np.asarray(dp)[:lengath]

IS_DRY = [0, 0, 0, 0, 0, 1, 0]
DOCONV_ALL = [0, 1, 1, 1, 1, 1, 1]
DOCONV_SKIP = [0, 1, 0, 1, 0, 1, 0]
TCONFIGS = {"t1": (0, DOCONV_ALL, 1800.0),
            "t2": (1, DOCONV_ALL, 3600.0),  # ztodt=2*dt; fixer fires
            "t3": (0, DOCONV_SKIP, 900.0)}

flux_args = (F(mu), F(md), F(du), F(eu), F(ed), F(dp))
out = dict(q=q, fracis=fracis, dpdry=dpdry, is_dry=np.array(IS_DRY),
           mu=np.asarray(mu), md=np.asarray(md), du=np.asarray(du),
           eu=np.asarray(eu), ed=np.asarray(ed), dp=np.asarray(dp),
           jt=np.asarray(jt_g), mx=np.asarray(mx_g),
           ideep=np.asarray(ideep), lengath=np.array(lengath))

gset = np.asarray(ideep)[:lengath] - 1
idle = sorted(set(range(ncol)) - set(gset.tolist()))
for cfg, (microp, doconv, dt) in TCONFIGS.items():
    dqdt = np.asarray(dtt(microp, doconv, IS_DRY, F(q), *flux_args,
                          jt_g, mx_g, ideep, lengath, F(fracis),
                          F(dpdry), dt))
    out[f"dqdt_{cfg}"] = dqdt
    out[f"doconvtran_{cfg}"] = np.array(doconv)
    out[f"dt_{cfg}"] = np.array(dt)
    assert np.isfinite(dqdt).all()
    assert np.all(dqdt[idle] == 0.0), "untriggered columns must no-op"
    assert np.all(dqdt[:, :, 0] == 0.0)
    for m in range(ncnst):
        if not doconv[m]:
            assert np.all(dqdt[:, :, m] == 0.0), (cfg, m)
    # uniform tracer no-op; column-mass conservation (gathered dp).
    # NOTE: the zm_microp fixer (t2) deliberately ADDS mass when a
    # column cannot absorb a negative tendency, so conservation is
    # asserted only for the microp=F configs.
    assert np.abs(dqdt[:, :, 4]).max() == 0.0
    if not microp:
        for m in (1, 2, 3, 5):
            if not doconv[m]:
                continue
            dpm = dpdry[:lengath] if IS_DRY[m] \
                else np.asarray(dp)[:lengath]
            col = np.abs(sum(dqdt[gset, k, m] * dpm[:, k]
                             for k in range(nlev)))
            scale = max(np.abs(dqdt[gset, :, m] * dpm).max(), 1e-30)
            assert col.max() < 1e-10 * scale, (cfg, m, col.max(), scale)
    print(f"cfg {cfg}: max|dqdt| per tracer =",
          [f"{np.abs(dqdt[:, :, m]).max():.2e}" for m in range(ncnst)])

# the fixer must actually fire and enforce positivity for m=7
dt2 = float(out["dt_t2"])
qnew = q[:, :, 6] + dt2 * out["dqdt_t2"][:, :, 6]
qraw = q[:, :, 6] + dt2 * out["dqdt_t1"][:, :, 6]
print(f"fixer: min q_new(m=7) without={qraw.min():.2e} "
      f"with={qnew.min():.2e}")
assert qraw.min() < 0.0, "m=7 must go negative without the fixer"
assert qnew.min() >= -1e-20

# --- momentum ---
winds1 = np.zeros((ncol, nlev, 2))
winds1[:, :, 0] = 25.0 * np.exp(-0.5 * ((sig - 0.25) / 0.18)**2) \
    + rng.uniform(-2.0, 2.0, (ncol, nlev))
winds1[:, :, 1] = 6.0 * np.sin(5.0 * sig) \
    + rng.uniform(-1.0, 1.0, (ncol, nlev))
winds2 = np.zeros((ncol, nlev, 2))
winds2[:, :, 0] = 7.0
winds2[:, :, 1] = -3.0
MNAMES = ("wind_tend", "pguall", "pgdall", "icwu", "icwd", "seten")
for cfg, w, dt in (("mom1", winds1, 1800.0), ("mom2", winds2, 1800.0)):
    res = dtm(lengath, F(w), *flux_args, jt_g, mx_g, ideep, dt)
    r = dict(zip(MNAMES, (np.asarray(v) for v in res)))
    out.update({f"{n}_{cfg}": v for n, v in r.items()})
    out[f"wind_in_{cfg}"] = w
    out[f"dt_{cfg}"] = np.array(dt)
    assert all(np.isfinite(v).all() for v in r.values())
    assert np.all(r["wind_tend"][idle] == 0.0)
    assert np.all(r["seten"][idle] == 0.0)
    assert np.array_equal(r["icwu"][idle], w[idle])
    for m in (0, 1):  # momentum column conservation
        col = np.abs(sum(r["wind_tend"][gset, k, m]
                         * np.asarray(dp)[:lengath, k]
                         for k in range(nlev)))
        assert col.max() < 1e-13, (cfg, m, col.max())
    print(f"cfg {cfg}: max|utend|={np.abs(r['wind_tend'][:, :, 0]).max():.2e}"
          f" max|seten|={np.abs(r['seten']).max():.2e}")
assert np.abs(out["wind_tend_mom2"]).max() < 1e-16  # uniform no-op

sha = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                     text=True, cwd="/work/E3SM").stdout.strip()
meta = {"scheme": "zm_transport zm_transport_tracer + "
                  "zm_transport_momentum (ZM convective transport)",
        "source_sha": sha, "ncol": ncol, "nlev": nlev,
        "ncnst": ncnst, "index_base": 1,
        "mass_flux_source": "eam_zm_conv_f drv_zm_conv_main re-run on "
                            "zm_conv_golden inputs, config a (EAMv3 "
                            "phys='default', is_first_step=T, "
                            "zm_microp=F); bitwise-equal to "
                            "zm_conv_golden *_a gathered fields",
        "zm_conv_params": {**p},
        "tracer_configs": {c: {"zm_microp": v[0], "doconvtran": v[1],
                               "dt": v[2]} for c, v in TCONFIGS.items()},
        "is_dry": IS_DRY,
        "momentum": {"momcu": 0.4, "momcd": 0.4, "nwind": 2,
                     "note": "compile-time parameters, no namelist "
                             "switch"},
        "constants": {"mbsth": 1e-15, "small": 1e-36,
                      "cdifr_min": 1e-6, "maxc_factor": 1e-12,
                      "flux_factor": 1e-12}}
gold = ROOT / "golden" / "zm_transport_golden.npz"
np.savez_compressed(gold, __metadata__=json.dumps(meta), **out)
print(f"wrote {gold}")
