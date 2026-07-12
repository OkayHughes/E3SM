"""Tier-1 golden replay + Tier-0 properties for the ZM deep-convection
main routine (eam_jax.zm_conv).

Golden archive: harness/gen_zm_conv_golden.py (run in the scream-dev
container) - 42 columns x 72 levels, 6 regime families, 3
configurations (a: EAMv3 defaults first step, b: EAMv3 defaults DCAPE
step, c: legacy flags), all with zm_microp=.false. (scope note in
PORTING_PLAN.md / eam_jax/zm_conv.py PORT_NOTES). Level/column
indices in the archive are 1-based (Fortran); the port is 0-based,
hence the -1 shifts.

Tolerances: everything replays at <= 1e-10 relative except `heat`
(observed max 6.0e-10, asserted at 2e-9). Documented reason: the
updraft profile mu = (exp(lambda*dz)-1)/dz and the entrainment rate
(rmue - mu)/dz are catastrophic cancellations of exp() outputs, and
the container glibc and host XLA exp() differ by 1 ulp on ~5% of
arguments (measured); the resulting few-ulp noise in mu/entr
propagates through the s_upd recursion into the dsdt flux difference
mu*(su - shat), which itself nearly cancels at levels of weak net
heating. qtnd is insensitive because its flux difference uses q_upd,
whose recursion has no exp-derived cancellation of the same size.
`ql`/`rprd` carry the same noise at denormal-scale entries, absorbed
by absolute floors far below physical significance.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
GOLDEN = ROOT / "golden" / "zm_conv_golden.npz"

from eam_jax import zm_conv  # noqa: E402

CFGS = ["a", "b", "c"]
NPC = 7        # columns per regime family (see generator)
GRAV = 9.80616
LATVAP = 2.501e6


@pytest.fixture(scope="module")
def gold():
    if not GOLDEN.exists():
        pytest.skip("golden archive missing (run harness in container)")
    return np.load(GOLDEN)


@pytest.fixture(scope="module")
def results(gold):
    """Run the port once per golden configuration."""
    meta = json.loads(str(gold["__metadata__"]))
    zc = zm_conv.make_zm_const()
    out = {}
    for cfg in CFGS:
        p = meta["params"][cfg]
        zp = zm_conv.make_zm_param(
            tau=p["tau"], alfa=p["alfa"], ke=p["ke"],
            dmpdz=p["dmpdz"], tpert_fix=bool(p["tpert_fix"]),
            tpert_fac=p["tpert_fac"], tiedke_add=p["tiedke_add"],
            c0_lnd=p["c0_lnd"], c0_ocn=p["c0_ocn"],
            num_cin=p["num_cin"], limcnv=p["limcnv"],
            mx_bot_lyr_adj=p["mx_bot_lyr_adj"],
            trig_dcape=bool(p["trig_dcape"]),
            trig_ull=bool(p["trig_ull"]),
            clos_dyn_adj=bool(p["clos_dyn_adj"]),
            no_deep_pbl=bool(p["no_deep_pbl"]),
            old_snow=bool(p["old_snow"]))
        out[cfg] = zm_conv.zm_conv_main(
            gold["t"], gold["q"], gold["omega"], gold["pmid"],
            gold["pint"], gold["pdel"], gold["geos"], gold["zmid"],
            gold["zint"], gold["pbl_hgt"], gold["tpert"],
            gold["landfrac"], gold["t_star"], gold["q_star"],
            float(gold["time_step"]), bool(p["is_first_step"]),
            zc, zp)
    return out


# ---------------------------------------------------------------------------
# Tier-1 golden replay
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("cfg", CFGS)
def test_gather_and_level_indices(gold, results, cfg):
    r = results[cfg]
    ng = r["lengath"]
    assert ng == int(gold[f"lengath_{cfg}"])
    np.testing.assert_array_equal(
        np.asarray(r["gather_index"][:ng]),
        gold[f"gather_index_{cfg}"][:ng] - 1)
    for f in ("msemax_klev_g", "jctop", "jcbot", "jt"):
        np.testing.assert_array_equal(np.asarray(r[f]),
                                      gold[f"{f}_{cfg}"] - 1)


# (field, rtol, atol); atol floors sit orders of magnitude below the
# physically meaningful range of each field (see module docstring)
REPLAY_SPECS = [
    ("cape", 1e-10, 1e-8),
    ("dcape", 1e-10, 1e-12),
    ("prec", 1e-10, 1e-18),
    ("rliq", 1e-10, 1e-18),
    ("qtnd", 1e-10, 1e-17),
    ("heat", 2e-9, 1e-13),     # cross-libm exp(); see docstring
    ("mcon", 1e-10, 1e-15),
    ("pflx", 1e-10, 1e-14),
    ("zdu", 1e-10, 1e-16),
    ("mflx_up", 1e-10, 1e-15),
    ("entr_up", 1e-10, 1e-16),
    ("detr_up", 1e-10, 1e-16),
    ("mflx_dn", 1e-10, 1e-15),
    ("entr_dn", 1e-10, 1e-16),
    ("p_del", 1e-13, 0.0),
    ("dsubcld", 1e-13, 0.0),
    ("ql", 1e-9, 1e-13),       # denormal-scale entries; see docstring
    ("rprd", 1e-9, 1e-16),
    ("dlf", 1e-10, 1e-16),
]


@pytest.mark.parametrize("cfg", CFGS)
@pytest.mark.parametrize("field,rtol,atol",
                         REPLAY_SPECS,
                         ids=[s[0] for s in REPLAY_SPECS])
def test_replay(gold, results, cfg, field, rtol, atol):
    np.testing.assert_allclose(np.asarray(results[cfg][field]),
                               gold[f"{field}_{cfg}"],
                               rtol=rtol, atol=atol)


# ---------------------------------------------------------------------------
# Tier-0 properties
# ---------------------------------------------------------------------------
def test_no_trigger_columns_are_exact_noops(results):
    """Columns outside the gather set (incl. the whole stable family)
    must produce exactly zero tendencies and precip."""
    for cfg in CFGS:
        r = results[cfg]
        ng = r["lengath"]
        active = set(np.asarray(r["gather_index"][:ng]).tolist())
        idle = sorted(set(range(len(np.asarray(r["prec"])))) - active)
        assert idle, "expected some inactive columns"
        for f in ("prec", "rliq"):
            np.testing.assert_array_equal(
                np.asarray(r[f])[idle], 0.0)
        for f in ("heat", "qtnd", "rprd", "dlf", "zdu"):
            np.testing.assert_array_equal(
                np.asarray(r[f])[idle], 0.0)
        # the stable-midlat family (cols 14-20) never triggers
        assert not (set(range(2 * NPC, 3 * NPC)) & active)


def test_precip_nonnegative(results):
    for cfg in CFGS:
        assert float(np.min(np.asarray(results[cfg]["prec"]))) >= 0.0
        assert float(np.min(np.asarray(results[cfg]["rliq"]))) >= 0.0


def test_water_balance(gold, results):
    """prec + rliq == -column integral of qtnd (to roundoff): the
    scheme removes exactly the water it dries the column by."""
    pdel = gold["pdel"]
    for cfg in CFGS:
        r = results[cfg]
        prec = np.asarray(r["prec"])
        rliq = np.asarray(r["rliq"])
        qint = (np.asarray(r["qtnd"]) * pdel).sum(axis=1) / GRAV / 1000.0
        assert np.abs(prec + rliq + qint).max() < 1e-18


def test_moist_enthalpy_balance(gold, results):
    """Column-integrated heating balances the latent heat of removed
    water: int(heat)dp/g == -L*int(qtnd)dp/g == L*rhow*(prec+rliq),
    to solver roundoff (~1e-15 relative in the Fortran itself)."""
    pdel = gold["pdel"]
    for cfg in CFGS:
        r = results[cfg]
        eint = (np.asarray(r["heat"]) * pdel).sum(axis=1) / GRAV
        lq = LATVAP * (np.asarray(r["qtnd"]) * pdel).sum(axis=1) / GRAV
        lp = LATVAP * 1000.0 * (np.asarray(r["prec"])
                                + np.asarray(r["rliq"]))
        scale = np.maximum(np.abs(eint), 1.0)
        assert (np.abs(eint + lq) / scale).max() < 1e-13
        assert (np.abs(eint - lp) / scale).max() < 1e-13


def test_mass_flux_consistency(results):
    """mcon (net convective mass flux) equals mflx_up + mflx_dn on the
    gathered rows; up-flux is nonnegative, down-flux nonpositive."""
    for cfg in CFGS:
        r = results[cfg]
        ng = r["lengath"]
        g = np.asarray(r["gather_index"][:ng])
        mu = np.asarray(r["mflx_up"])[:ng]
        md = np.asarray(r["mflx_dn"])[:ng]
        net = np.asarray(r["mcon"])[g, :mu.shape[1]]
        # (a+b)*cbmf vs a*cbmf + b*cbmf: roundoff-level only
        np.testing.assert_allclose(net, mu + md, rtol=1e-13,
                                   atol=1e-20)
        assert mu.min() >= 0.0
        assert md.max() <= 0.0


def test_plume_confined_between_top_and_base(results):
    """Updraft mass flux vanishes at/above the plume top and below the
    launch level."""
    for cfg in CFGS:
        r = results[cfg]
        ng = r["lengath"]
        mu = np.asarray(r["mflx_up"])[:ng]
        jt = np.asarray(r["jt"])[:ng]
        jb = np.asarray(r["msemax_klev_g"])[:ng]
        pver = mu.shape[1]
        for i in range(ng):
            assert (mu[i, :jt[i] + 1] == 0.0).all()
            assert (mu[i, jb[i] + 1:pver] == 0.0).all()
