"""Tier-1.5 chained replay of the full ZM deep-convection sequence
(eam_jax.zm_intr.zm_tend vs the Fortran chain).

Golden archive: harness/gen_zm_chain_golden.py (scream-dev container)
-- 42 columns x 72 levels, 3 consecutive pseudo-timesteps (dt=1800 s)
of the complete zm_conv_tend sequence (DCAPE-triggered zm_conv_main ->
physics_update -> zm_conv_evap -> physics_update -> momentum + 2-tracer
transport -> summed ptend), with tendency feedback into the state, a
fixed destabilizing inter-step forcing (dynamics stand-in), end-of-step
T_STAR/Q_STAR recording for the DCAPE trigger, and a geopotential_t
refresh of zmid/zint. The JAX side replays the identical chain from the
stored initial state through zm_tend + eam_jax.geopotential and is
compared per step (trigger set exact, tendencies at diagnostic
tolerances) and on the final state.

Final-state tolerances (target 1e-10 relative; measured on this
archive, container libm vs host XLA):
  t     1.6e-12   zmid  7.2e-14   tr2  4.6e-11    -> asserted 1e-10
  qv    1.1e-10   tr1   1.9e-10   u    2.2e-10   v 1.6e-10
                                                  -> asserted 1e-9
The four looser fields inherit the documented zm_conv `heat` exception
(tests/test_zm_conv.py: 1-ulp cross-libm exp() through the
(exp(x)-1)/x mass-flux cancellation, bounded there at 2e-9 per call):
per-step s_tend noise is measured at 7.6e-10/3.4e-9/2.2e-9 over the
three steps, and each step's noise perturbs the next step's trigger
inputs (t via s_tend/cpair, the winds via the mu-dependent momentum
transport, qv/tr1 via the plume fluxes), so the accumulated libm ulp
difference through 3 chained steps lands at ~2e-10 -- feedback
amplification of documented kernel-level noise, not an interface
discrepancy (the gather sets, level indices, prec and snow stay
exact/near-exact throughout: prec <= 1.9e-12 rel, snow bitwise 0).
"""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
GOLDEN = ROOT / "golden" / "zm_chain_golden.npz"

from eam_jax import zm_conv, zm_intr  # noqa: E402
from eam_jax.constants import CPAIR, GRAVIT, RAIR, ZVIR  # noqa: E402
from eam_jax.geopotential import geopotential_t  # noqa: E402


@pytest.fixture(scope="module")
def gold():
    if not GOLDEN.exists():
        pytest.skip("golden archive missing (run harness in container)")
    return np.load(GOLDEN)


@pytest.fixture(scope="module")
def chain(gold):
    """Run the 3-step JAX chain; returns (per-step results, final
    state dict)."""
    meta = json.loads(str(gold["__metadata__"]))
    P = meta["params"]
    dt = meta["dt"]
    nsteps = meta["nsteps"]
    zc = zm_conv.make_zm_const()
    zp = zm_conv.make_zm_param(
        tau=P["tau"], alfa=P["alfa"], ke=P["ke"], dmpdz=P["dmpdz"],
        tpert_fix=bool(P["tpert_fix"]), tpert_fac=P["tpert_fac"],
        tiedke_add=P["tiedke_add"], c0_lnd=P["c0_lnd"],
        c0_ocn=P["c0_ocn"], num_cin=P["num_cin"], limcnv=P["limcnv"],
        mx_bot_lyr_adj=P["mx_bot_lyr_adj"],
        trig_dcape=bool(P["trig_dcape"]),
        trig_ull=bool(P["trig_ull"]),
        clos_dyn_adj=bool(P["clos_dyn_adj"]),
        no_deep_pbl=bool(P["no_deep_pbl"]),
        old_snow=bool(P["old_snow"]))

    t = np.array(gold["t0"])
    qv = np.array(gold["qv0"])
    tr1 = np.array(gold["tr1_0"])
    tr2 = np.array(gold["tr2_0"])
    u = np.array(gold["u0"])
    v = np.array(gold["v0"])
    zmid = np.array(gold["zmid0"])
    zint = np.array(gold["zint0"])
    pmid, pint, pdel = gold["pmid"], gold["pint"], gold["pdel"]
    piln = np.log(pint)
    pmln = np.log(pmid)
    rpdel = 1.0 / pdel
    t_star, q_star = t.copy(), qv.copy()

    steps = []
    for step in range(nsteps):
        q3 = np.stack([qv, tr1, tr2], axis=2)
        r = zm_intr.zm_tend(
            t, q3, u, v, gold["omega"], pmid, pint, pdel,
            gold["geos"], zmid, zint, gold["pblh"], gold["tpert"],
            gold["landfrac"], gold["cld"], gold["fracis"],
            gold["doconvtran"].astype(bool), t_star, q_star,
            dt, step == 0, zc, zp)
        steps.append(r)

        # apply the summed ptend in one update (physpkg ptend_all),
        # with the qneg3-style clips of the golden update rule
        t = t + dt * np.asarray(r["s_tend"]) / CPAIR
        qt = np.asarray(r["q_tend"])
        qv = np.maximum(qv + dt * qt[:, :, 0], zm_intr.QMIN_VAPOR)
        tr1 = np.maximum(tr1 + dt * qt[:, :, 1], 0.0)
        tr2 = np.maximum(tr2 + dt * qt[:, :, 2], 0.0)
        u = u + dt * np.asarray(r["u_tend"])
        v = v + dt * np.asarray(r["v_tend"])

        # end-of-step DCAPE state, then forcing + geopotential refresh
        t_star, q_star = t.copy(), qv.copy()
        t = t + dt * gold["t_forc"]
        qv = qv + dt * gold["q_forc"]
        zi_new, zm_new = geopotential_t(piln, pmln, pint, pmid, pdel,
                                        rpdel, t, qv, RAIR, GRAVIT,
                                        ZVIR)
        zint = np.asarray(zi_new)
        zmid = np.asarray(zm_new)

    final = dict(t=t, qv=qv, tr1=tr1, tr2=tr2, u=u, v=v,
                 zmid=zmid, zint=zint, t_star=t_star, q_star=q_star)
    return steps, final


def _rel(a, b, floor):
    a, b = np.asarray(a), np.asarray(b)
    return float((np.abs(a - b) / np.maximum(np.abs(b), floor)).max())


# ---------------------------------------------------------------------------
# per-step interface checks
# ---------------------------------------------------------------------------
def test_trigger_set_exact(gold, chain):
    """The gather (trigger) set must match the Fortran exactly on all
    steps: any interface misunderstanding (t_star bookkeeping, DCAPE
    sign, state feedback) would change which columns convect."""
    steps, _ = chain
    for step, r in enumerate(steps):
        ng = r["lengath"]
        assert ng == int(gold[f"lengath_{step}"]), step
        np.testing.assert_array_equal(
            np.asarray(r["ideep"][:ng]),
            gold[f"gather_index_{step}"][:ng] - 1)


# (field-in-result, golden key stem, rel tol, floor); tolerances are
# ~3-5x the measured maxima (see module docstring)
STEP_SPECS = [
    ("s_tend", "s_tend", 1e-8, 1e-4),
    ("qv_tend", "qv_tend", 5e-9, 1e-12),
    ("u_tend", "u_tend", 5e-9, 1e-9),
    ("v_tend", "v_tend", 5e-9, 1e-9),
    ("tr1_tend", "dqdt_tr1", 5e-9, 1e-13),
    ("tr2_tend", "dqdt_tr2", 5e-9, 1e-13),
    ("prec", "prec", 1e-11, 1e-12),
    ("snow", "snow", 1e-11, 1e-12),
    ("cape", "cape", 1e-9, 1e-6),
    ("dcape", "dcape", 1e-8, 1e-8),
]


@pytest.mark.parametrize("field,gstem,rtol,floor", STEP_SPECS,
                         ids=[s[0] for s in STEP_SPECS])
def test_per_step_tendencies(gold, chain, field, gstem, rtol, floor):
    steps, _ = chain
    for step, r in enumerate(steps):
        if field == "qv_tend":
            a = np.asarray(r["q_tend"])[:, :, 0]
        elif field == "tr1_tend":
            a = np.asarray(r["q_tend"])[:, :, 1]
        elif field == "tr2_tend":
            a = np.asarray(r["q_tend"])[:, :, 2]
        else:
            a = np.asarray(r[field])
        d = _rel(a, gold[f"{gstem}_{step}"], floor)
        assert d < rtol, (field, step, d)


# ---------------------------------------------------------------------------
# final-state agreement
# ---------------------------------------------------------------------------
FINAL_SPECS = [
    # (field, rel tol, floor) -- measured maxima in the docstring
    ("t", 1e-10, 1.0),
    ("zmid", 1e-10, 1.0),
    ("zint", 1e-10, 1.0),
    ("tr2", 1e-10, 1e-12),
    ("qv", 1e-9, 1e-9),    # measured 1.1e-10; heat-noise feedback
    ("tr1", 1e-9, 1e-12),  # measured 1.9e-10
    ("u", 1e-9, 1e-3),     # measured 2.2e-10
    ("v", 1e-9, 1e-3),     # measured 1.6e-10
    ("t_star", 1e-10, 1.0),
    ("q_star", 1e-9, 1e-9),
]


@pytest.mark.parametrize("field,rtol,floor", FINAL_SPECS,
                         ids=[s[0] for s in FINAL_SPECS])
def test_final_state(gold, chain, field, rtol, floor):
    _, final = chain
    gkey = {"t_star": "tstar_final", "q_star": "qstar_final"}.get(
        field, f"{field}_final")
    d = _rel(final[field], gold[gkey], floor)
    assert d < rtol, (field, d)


def test_chain_stays_convecting(gold, chain):
    """The DCAPE trigger must keep the chain alive on every step (a
    dead chain would trivially 'agree')."""
    steps, _ = chain
    for r in steps:
        assert r["lengath"] > 0


def test_chain_water_accounting(gold, chain):
    """Per step: column water removed from vapor+tracers equals
    precip + reserved liquid (the zm_conv balance) minus what
    evaporation put back -- i.e. prec+rliq == -int(qv_tend)dp/g
    (transport conserves tracer mass separately)."""
    steps, _ = chain
    pdel = np.asarray(gold["pdel"])
    for step, r in enumerate(steps):
        prec = np.asarray(r["prec"])
        rliq = np.asarray(r["rliq"])
        qint = (np.asarray(r["q_tend"])[:, :, 0] * pdel).sum(1) \
            / GRAVIT / 1000.0
        assert np.abs(prec + rliq + qint).max() < 1e-17, step
