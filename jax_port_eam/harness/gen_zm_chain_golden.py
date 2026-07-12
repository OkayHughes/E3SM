#!/usr/bin/env python3
"""Tier-1.5 chained-replay golden for the full ZM deep-convection
sequence (run in the scream-dev container).

Drives the FULL zm_conv_tend chain through the f2py extensions for 3
consecutive pseudo-timesteps (dt = 1800 s) on the 42 zm_conv golden
columns, feeding tendencies back into the state exactly as the model
does, and capturing the state before/after every step. The JAX side
(tests/test_zm_chain.py) replays the same chain through
eam_jax.zm_intr.zm_tend and must reproduce the final state -- this
catches interface misunderstandings between the per-kernel ports
(gathered-vs-scattered layouts, mb/Pa conventions, update ordering,
DCAPE t_star/q_star bookkeeping) that kernel goldens cannot.

Per step, exactly the zm_conv_tend sequence (see eam_jax/zm_intr.py
PORT_NOTES): zm_conv_main (eam_zm_conv_f) -> physics_update (t1 =
t+dt*heat/cpair, q1 = max(q+dt*qtnd, qmin)) -> zm_conv_evap
(eam_zm_evap_f, prdprec=rprd, pbuf cloud fraction, prec inout) ->
physics_update (q2) -> zm_transport_momentum + zm_transport_tracer
(eam_zm_transport_f, 2 tracers, vapor never transported) -> summed
ptend (s = heat+evap_s+seten; qv = qtnd+evap_q) applied to the state
in ONE update, as physpkg applies ptend_all.

Between steps (standing in for dynamics + the rest of physics):
  - T_STAR/Q_STAR are recorded from the end-of-step state (physpkg
    records them at the end of tphysac; on step 0 zm_conv_tend seeds
    them with the current state, reproduced here);
  - a fixed destabilizing forcing (boundary-layer warming/moistening
    + mid-level cooling, stored in the archive) is applied so the
    DCAPE trigger (cape>0 AND dcape>0) stays alive after step 0;
  - zmid/zint are recomputed from the new (t, qv) with the real
    geopotential_t (eam_geopotential_f, SE branch), as physics_update
    does after every T/q change.

Scope: zm_microp=.false. (PORTING_PLAN row 9), old_snow=.true. (the
zm_param_t default -- the configuration consistent with microp off),
MCSP off (mcsp_enabled=F; EAMv3 phys="default" would set
zmconv_MCSP_heat_coeff=0.3 but zm_conv_mcsp is not ported), no
aerosol transport. Params otherwise EAMv3 phys="default"
(zmconv_ke=2.5e-6 dyn=se microphys=p3). qneg3 clip at qmin=1e-12 for
vapor and max(.,0) for tracers is applied on both sides (asserted to
never actually fire on this chain).
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "fmod"))
import eam_geopotential_f  # noqa: E402
import eam_zm_conv_f  # noqa: E402
import eam_zm_evap_f  # noqa: E402
import eam_zm_transport_f  # noqa: E402

F = lambda a: np.asarray(a, dtype=np.float64, order="F")  # noqa: E731

# constants (share/util/shr_const_mod.F90 / physconst)
GRAV = 9.80616
CPAIR = 1.00464e3
RGAS = 6.02214e26 * 1.38065e-23
RDAIR = RGAS / 28.966
RWV = RGAS / 18.016
ZVIR = RWV / RDAIR - 1.0        # physconst zvir (NOT the ZM 1.608)
QMIN_VAPOR = 1.0e-12
KE_V3 = 2.5e-6
DT = 1800.0
NSTEPS = 3
NCNST = 3                        # vapor + 2 transported tracers

dc = eam_zm_conv_f.zm_conv_driver
de = eam_zm_evap_f.zm_evap_driver
dtr = eam_zm_transport_f.zm_transport_driver
dg = eam_geopotential_f.geopotential_driver
for f in (dc.drv_zm_conv_main, de.drv_zm_conv_evap,
          dtr.drv_transport_tracer, dtr.drv_transport_momentum,
          dg.drv_geopotential_t):
    assert "array('d')" in f.__doc__
dc.drv_init()
de.drv_init()

# --- initial state: the zm_conv golden profiles ---
g = np.load(ROOT / "golden" / "zm_conv_golden.npz")
zmeta = json.loads(str(g["__metadata__"]))
ncol, nlev = np.asarray(g["t"]).shape
limcnv = int(g["limcnv"])
P = dict(tau=3600.0, alfa=0.14, ke=KE_V3, dmpdz=-0.7e-3,
         tpert_fix=1, tpert_fac=2.0, tiedke_add=0.8,
         c0_lnd=0.0020, c0_ocn=0.0020, num_cin=1,
         mx_bot_lyr_adj=1, trig_dcape=1, trig_ull=1,
         clos_dyn_adj=1, old_snow=1, no_deep_pbl=0)

pmid = np.asarray(g["pmid"]); pint = np.asarray(g["pint"])
pdel = np.asarray(g["pdel"]); geos = np.asarray(g["geos"])
omega = np.asarray(g["omega"])
pblh = np.asarray(g["pbl_hgt"]); tpert = np.asarray(g["tpert"])
landfrac = np.asarray(g["landfrac"])
piln = np.log(pint); pmln = np.log(pmid); rpdel = 1.0 / pdel
sig = pmid / pmid[:, -1:]

rng = np.random.default_rng(20260715)
t = np.asarray(g["t"]).copy()
qv = np.asarray(g["q"]).copy()
u = 20.0 * np.exp(-0.5 * ((sig - 0.3) / 0.2) ** 2) \
    + rng.uniform(-1.5, 1.5, (ncol, nlev))
v = 5.0 * np.sin(4.0 * sig) + rng.uniform(-1.0, 1.0, (ncol, nlev))
tr1 = 1.0e-6 * np.exp(4.0 * (sig - 1.0)) \
    * (1.0 + 0.3 * rng.uniform(-1, 1, (ncol, nlev)))       # BL tracer
tr2 = 1.0e-9 + 2.0e-6 * np.exp(-0.5 * ((sig - 0.45) / 0.12) ** 2)
cld = np.clip(0.25 + 0.5 * np.exp(-0.5 * ((sig - 0.55) / 0.2) ** 2)
              * (1.0 + 0.3 * rng.uniform(-1, 1, (ncol, nlev))),
              0.0, 0.9)
fracis = np.ones((ncol, nlev, NCNST))
fracis[:, :, 2] = 0.3 + 0.7 * sig
DOCONV = np.array([0, 1, 1])
IS_DRY = np.array([0, 0, 0])
zmid = np.asarray(g["zmid"]).copy()
zint = np.asarray(g["zint"]).copy()

# fixed destabilizing forcing (dynamics stand-in), from the initial
# state: ~0.36 K/step BL warming, mid-level cooling, 1.5%/step BL
# moistening
t_forc = (2.0e-4 * np.clip((sig - 0.8) / 0.2, 0.0, 1.0)
          - 4.0e-5 * np.exp(-0.5 * ((sig - 0.5) / 0.15) ** 2))
q_forc = 0.015 / 1800.0 * qv * np.clip((sig - 0.75) / 0.25, 0.0, 1.0)

CNAMES = ("lengath", "gather_index", "msemax_klev_g", "jctop", "jcbot",
          "jt", "prec", "heat", "qtnd", "cape", "dcape", "mcon", "pflx",
          "zdu", "mflx_up", "entr_up", "detr_up", "mflx_dn", "entr_dn",
          "p_del", "dsubcld", "ql", "rliq", "rprd", "dlf")
ENAMES = ("tend_s", "tend_q", "tend_s_snwprd", "tend_s_snwevmlt",
          "prec", "snow", "ntprprd", "ntsnprd", "flxprec", "flxsnow")
MNAMES = ("wind_tend", "pguall", "pgdall", "icwu", "icwd", "seten")

out = dict(t0=t.copy(), qv0=qv.copy(), tr1_0=tr1.copy(),
           tr2_0=tr2.copy(), u0=u.copy(), v0=v.copy(),
           zmid0=zmid.copy(), zint0=zint.copy(),
           omega=omega, pmid=pmid, pint=pint, pdel=pdel, geos=geos,
           pblh=pblh, tpert=tpert, landfrac=landfrac, cld=cld,
           fracis=fracis, doconvtran=DOCONV, is_dry=IS_DRY,
           t_forc=t_forc, q_forc=q_forc,
           dt=np.array(DT), nsteps=np.array(NSTEPS),
           limcnv=np.array(limcnv))

t_star = t.copy()
q_star = qv.copy()
zeros2 = np.zeros((ncol, nlev))
for step in range(NSTEPS):
    is_first = 1 if step == 0 else 0
    out[f"tstar_{step}"] = t_star.copy()
    out[f"qstar_{step}"] = q_star.copy()
    out[f"t_in_{step}"] = t.copy()
    out[f"qv_in_{step}"] = qv.copy()
    out[f"tr1_in_{step}"] = tr1.copy()
    out[f"tr2_in_{step}"] = tr2.copy()
    out[f"u_in_{step}"] = u.copy()
    out[f"v_in_{step}"] = v.copy()
    out[f"zmid_in_{step}"] = zmid.copy()
    out[f"zint_in_{step}"] = zint.copy()

    # 1. zm_conv_main
    res = dc.drv_zm_conv_main(
        DT, is_first, limcnv, P["no_deep_pbl"], P["tau"], P["alfa"],
        P["ke"], P["dmpdz"], P["tpert_fix"], P["tpert_fac"],
        P["tiedke_add"], P["c0_lnd"], P["c0_ocn"], P["num_cin"],
        P["mx_bot_lyr_adj"], P["trig_dcape"], P["trig_ull"],
        P["clos_dyn_adj"], P["old_snow"],
        F(t), F(qv), F(omega), F(pmid), F(pint), F(pdel), geos,
        F(zmid), F(zint), pblh, tpert, landfrac, F(t_star), F(q_star))
    c = dict(zip(CNAMES, (np.asarray(x) for x in res)))
    lengath = int(c["lengath"])
    assert lengath > 0, f"step {step}: chain died (no convection)"

    # 3. physics_update for 'zm_conv_main' (qneg3 must be inert)
    t1 = t + DT * c["heat"] / CPAIR
    q1_raw = qv + DT * c["qtnd"]
    assert q1_raw.min() >= QMIN_VAPOR, "qneg3 fired: redesign forcing"
    q1 = np.maximum(q1_raw, QMIN_VAPOR)

    # 4. zm_conv_evap on state1
    er = dict(zip(ENAMES, (np.asarray(x) for x in de.drv_zm_conv_evap(
        KE_V3, P["old_snow"], DT, F(pmid), F(pdel), F(t1), F(q1),
        F(c["rprd"]), F(cld), F(zeros2), F(zeros2), c["prec"]))))

    # 5. physics_update for 'zm_conv_evap'
    q2_raw = q1 + DT * er["tend_q"]
    assert q2_raw.min() >= QMIN_VAPOR
    q2 = np.maximum(q2_raw, QMIN_VAPOR)

    # 6. momentum transport (u/v untouched so far)
    winds = np.stack([u, v], axis=2)
    mres = dtr.drv_transport_momentum(
        lengath, F(winds), F(c["mflx_up"]), F(c["mflx_dn"]),
        F(c["detr_up"]), F(c["entr_up"]), F(c["entr_dn"]),
        F(c["p_del"]), c["jt"], c["msemax_klev_g"], c["gather_index"],
        DT)
    m = dict(zip(MNAMES, (np.asarray(x) for x in mres)))

    # 7. tracer transport on state1%q (vapor slot updated, never
    # transported)
    q3 = np.stack([q2, tr1, tr2], axis=2)
    dqdt = np.asarray(dtr.drv_transport_tracer(
        0, DOCONV, IS_DRY, F(q3), F(c["mflx_up"]), F(c["mflx_dn"]),
        F(c["detr_up"]), F(c["entr_up"]), F(c["entr_dn"]),
        F(c["p_del"]), c["jt"], c["msemax_klev_g"], c["gather_index"],
        lengath, F(fracis), F(zeros2), DT))

    # 8. summed ptend, applied in one update (as physpkg applies
    # ptend_all); order of the sum matches physics_ptend_sum calls
    s_tend = (c["heat"] + er["tend_s"]) + m["seten"]
    qv_tend = c["qtnd"] + er["tend_q"]
    out[f"s_tend_{step}"] = s_tend
    out[f"qv_tend_{step}"] = qv_tend
    out[f"u_tend_{step}"] = m["wind_tend"][:, :, 0]
    out[f"v_tend_{step}"] = m["wind_tend"][:, :, 1]
    out[f"dqdt_tr1_{step}"] = dqdt[:, :, 1]
    out[f"dqdt_tr2_{step}"] = dqdt[:, :, 2]
    out[f"lengath_{step}"] = np.array(lengath)
    out[f"gather_index_{step}"] = c["gather_index"]
    out[f"cape_{step}"] = c["cape"]
    out[f"dcape_{step}"] = c["dcape"]
    out[f"prec_{step}"] = er["prec"]
    out[f"snow_{step}"] = er["snow"]
    out[f"rliq_{step}"] = c["rliq"]
    out[f"heat_{step}"] = c["heat"]
    out[f"qtnd_{step}"] = c["qtnd"]
    out[f"evap_s_{step}"] = er["tend_s"]
    out[f"evap_q_{step}"] = er["tend_q"]
    out[f"seten_{step}"] = m["seten"]
    out[f"flxprec_{step}"] = er["flxprec"]
    out[f"flxsnow_{step}"] = er["flxsnow"]

    t = t + DT * s_tend / CPAIR
    qv_raw = qv + DT * qv_tend
    tr1_raw = tr1 + DT * dqdt[:, :, 1]
    tr2_raw = tr2 + DT * dqdt[:, :, 2]
    assert qv_raw.min() >= QMIN_VAPOR, "vapor qneg3 fired"
    # the qneg3-style tracer clip at 0 IS allowed to fire (it does,
    # on a single point: the min(chat,const) flux limiter bounds the
    # tendency but does not guarantee positivity over a full ztodt);
    # both chain sides apply the identical clip
    nclip = int((tr1_raw < 0.0).sum() + (tr2_raw < 0.0).sum())
    print(f"  step {step}: tracer clip points = {nclip}")
    qv = np.maximum(qv_raw, QMIN_VAPOR)
    tr1 = np.maximum(tr1_raw, 0.0)
    tr2 = np.maximum(tr2_raw, 0.0)
    u = u + DT * m["wind_tend"][:, :, 0]
    v = v + DT * m["wind_tend"][:, :, 1]

    out[f"t_phys_{step}"] = t.copy()
    out[f"qv_phys_{step}"] = qv.copy()

    # end-of-step T_STAR/Q_STAR (physpkg, end of tphysac)
    t_star = t.copy()
    q_star = qv.copy()

    # dynamics stand-in forcing + geopotential refresh
    t = t + DT * t_forc
    qv = qv + DT * q_forc
    rair2 = np.full((ncol, nlev), RDAIR)
    zvir2 = np.full((ncol, nlev), ZVIR)
    zi_new, zm_new = dg.drv_geopotential_t(
        F(piln), F(pmln), F(pint), F(pmid), F(pdel), F(rpdel),
        F(t), F(qv), F(rair2), GRAV, F(zvir2))
    zint = np.asarray(zi_new)
    zmid = np.asarray(zm_new)

    print(f"step {step}: lengath={lengath} "
          f"prec_max[mm/day]={er['prec'].max() * 8.64e7:.2f} "
          f"|dT|max={np.abs(DT * s_tend / CPAIR).max():.3f} K "
          f"dcape_max={c['dcape'].max():.3f}")

out.update(t_final=t, qv_final=qv, tr1_final=tr1, tr2_final=tr2,
           u_final=u, v_final=v, zmid_final=zmid, zint_final=zint,
           tstar_final=t_star, qstar_final=q_star)

sha = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                     text=True, cwd="/work/E3SM").stdout.strip()
meta = {"scheme": "Tier-1.5 chained replay: zm_conv_main -> "
                  "physics_update -> zm_conv_evap -> physics_update "
                  "-> zm_transport_momentum -> zm_transport_tracer, "
                  "3 steps with tendency feedback",
        "source_sha": sha, "ncol": ncol, "nlev": nlev,
        "ncnst": NCNST, "nsteps": NSTEPS, "dt": DT,
        "params": {**P, "limcnv": limcnv, "zm_microp": 0,
                   "mcsp_enabled": 0},
        "old_snow_note": "old_snow=T is the zm_param_t default and "
                         "the configuration consistent with "
                         "zm_microp=F (zm_conv_intr only sets it F "
                         "under zmconv_microp)",
        "ke_note": "zmconv_ke=2.5e-6 (namelist_defaults_eam.xml "
                   "phys='default' dyn='se' microphys='p3')",
        "update_rule": "state += dt*sum(ptend) in one update "
                       "(physpkg ptend_all); vapor qneg3 clip at "
                       "1e-12 and tracer clip at 0 asserted inert",
        "interstep": "t_star/q_star recorded at end of step "
                     "(tphysac); then fixed destabilizing forcing "
                     "(t_forc, q_forc arrays); then geopotential_t "
                     "(SE branch) refresh of zmid/zint with "
                     "physconst zvir",
        "constants": {"cpair": CPAIR, "grav": GRAV, "rdair": RDAIR,
                      "zvir_geopotential": ZVIR,
                      "qmin_vapor": QMIN_VAPOR}}
gold = ROOT / "golden" / "zm_chain_golden.npz"
np.savez_compressed(gold, __metadata__=json.dumps(meta), **out)
print(f"wrote {gold}")
