"""Golden replay + Tier-0 property tests for the ZM convective
microphysics port (eam_jax/zm_microphysics.py) and the zm_microp=True
extension of zm_conv_main (eam_jax/zm_conv.py).

Golden: harness/gen_zm_microp_golden.py run in the scream-dev
container against the REAL zm_microphysics.F90 + activate_drop_mam +
nucleate_ice_conv (modal aerosols; see PORTING_PLAN.md row 9).

Tolerances: target 1e-12 relative. Measured exceptions are documented
at the assertions: the kernel chains glibc-vs-XLA/cephes libm calls
(pow with real exponents, exp, erf, tgamma) through iterated plume
recursions and conservation-ratio branches, which amplifies 1-ulp
libm differences; absolute floors sit at the denormal residue of
cancelling process-rate sums.
"""
import json
from pathlib import Path

import numpy as np
import pytest

from eam_jax import zm_microphysics as zmp
from eam_jax.zm_cape import make_zm_const

GOLD = Path(__file__).resolve().parents[1] / "golden" / "zm_microp_golden.npz"
D = np.load(GOLD)
META = json.loads(str(D["__metadata__"]))
STATE_FIELDS = META["state_fields"]
DIAG_FIELDS = META["diag_fields"]
MICROP_FIELDS = META["microp_fields"]
P = META["params"]
AERO_META = META["aero"]


def _aero_dict(numg, mmrg, dgnumg):
    m_acc, m_ait, m_crs = AERO_META["mode_idx"]
    l_dst, l_ncl, l_so4 = AERO_META["coarse_species_idx"]
    return dict(scheme="modal",
                nspec=AERO_META["nspec"],
                mode_accum=m_acc - 1, mode_aitken=m_ait - 1,
                mode_coarse=m_crs - 1, coarse_dust=l_dst - 1,
                coarse_nacl=l_ncl - 1, coarse_so4=l_so4 - 1,
                sigmag_aitken=AERO_META["sigmag"][m_ait - 1],
                sigmag_amode=np.array(AERO_META["sigmag"]),
                specdens=D["aero_specdens"],
                spechygro=D["aero_spechygro"],
                voltonumblo=D["aero_voltonumblo"],
                voltonumbhi=D["aero_voltonumbhi"],
                numg=numg, mmrg=mmrg, dgnumg=dgnumg)


ZC = make_zm_const()


@pytest.fixture(scope="module")
def kernel_result():
    aero = _aero_dict(D["k_numg"], D["k_mmrg"], D["k_dgnum"])
    out = zmp.zm_mphy(
        msg=int(D["k_msg"]),
        jb=D["k_jb"] - 1, jt=D["k_jt"] - 1, jlcl=D["k_jlcl"] - 1,
        su=D["k_su"], qu=D["k_qu"], mu=D["k_mu"], du=D["k_du"],
        eu=D["k_eu"], zf=D["k_zf"], pm=D["k_pm"], te=D["k_te"],
        qe=D["k_qe"], gamhat=D["k_gamhat"], eps0=D["k_eps0"],
        cmel=D["k_cmel"], cmei=D["k_cmei"], aero=aero,
        deltat=P["deltat"], auto_fac=P["auto_fac"],
        accr_fac=P["accr_fac"], dcs=P["micro_dcs"],
        grav=ZC["grav"], cp=ZC["cpair"], rd=ZC["rdair"],
        lamc0=P["lambdadpcu0"], pgam0=P["mudpcu0"])
    return {k: np.asarray(v) for k, v in out.items()}


def _cmp(name, got, want, rtol, atol):
    err = np.abs(got - want)
    tol = atol + rtol * np.abs(want)
    bad = err > tol
    if np.any(bad):
        i = np.unravel_index(np.argmax(err / np.maximum(tol, 1e-300)),
                             err.shape)
        raise AssertionError(
            f"{name}: {bad.sum()} pts out of tol; worst at {i}: "
            f"got {got[i]!r} want {want[i]!r} "
            f"abs {err[i]:.3e} rel {err[i] / max(abs(want[i]), 1e-300):.3e}")


# measured replay agreement, config k (see module docstring).
# Default 1e-12: qc/qi/nc/ni/wu/lamc/pgam/frz + most process diags
# replay below it. Documented exceptions (all measured, 2026-07):
# - precip species chain (qr/nr/qni/ns/qg/ng + their detrainment and
#   diag partners): sub-QSMALL cancellation residues in the vertical-
#   integration sums (mu*q + dz*(tend+fall)*arcf ~ 1e-24) are
#   amplified through fallspeed pow() chains and conservation-ratio
#   branches; measured <= 2.1e-10 rel -> rtol 5e-10.
# - fhmrm: homogeneous rain freezing at T=233.15K acts on qr values
#   that are pure cancellation residue (~1e-17); their SIGN decides
#   the branch, so one side can freeze a denormal qr the other
#   clipped to zero. Measured abs <= 1.3e-10 with ZERO measured
#   impact on frz/sprd/qg (<= 1e-21 abs) -> atol 5e-10.
KTOL = {"default": (1e-12, 0.0)}
for _nm in ("qr", "nr", "qni", "ns", "qg", "ng", "qnide", "nsde",
            "sprd", "rprd", "accrlm", "accrln", "fallrm", "fallrn",
            "fallsm", "fallsn", "fallgm", "fallgn", "dsfm", "dsfn",
            "accgrm", "accgrn", "accgsrm", "accgsrn", "fhtimm",
            "fhtctm", "autolm", "autoln", "accslm", "accsln"):
    KTOL[_nm] = (5e-10, 1e-22)
KTOL["fhmrm"] = (5e-10, 5e-10)
# bergnn (Bergeron nprb ~ prb*ncic/qcic libm chain) and trspcm (the
# mu*qc(k)-mu*qc(k-1) flux difference cancels to ~1e-22) measured at
# 2.2e-12 / 1.8e-12 rel
KTOL["bergnn"] = (5e-12, 0.0)
KTOL["trspcm"] = (5e-12, 1e-21)


class TestKernelReplay:
    @pytest.mark.parametrize("i", range(len(STATE_FIELDS)))
    def test_state(self, kernel_result, i):
        nm = STATE_FIELDS[i]
        want = D["k_state"][:, :, i]
        rtol, atol = KTOL.get(nm, KTOL["default"])
        _cmp(nm, kernel_result[nm], want, rtol, atol)

    @pytest.mark.parametrize("i", range(len(DIAG_FIELDS)))
    def test_diag(self, kernel_result, i):
        nm = DIAG_FIELDS[i]
        want = D["k_diag"][:, :, i]
        rtol, atol = KTOL.get(nm, KTOL["default"])
        _cmp(nm, kernel_result[nm], want, rtol, atol)


# ---------------------------------------------------------------------------
# Tier-1 golden replay: full zm_conv_main with zm_microp=True
# ---------------------------------------------------------------------------
from eam_jax import zm_conv  # noqa: E402

MAIN_CFGS = ["m", "n"]
# golden zm_microp_st field name -> port key in out["microp"]
MICROP_NAME_MAP = {"qliq": "qc", "qice": "qi", "qrain": "qr",
                   "qsnow": "qni", "qgraupel": "qg", "qnl": "nc",
                   "qni": "ni", "qnr": "nr", "qns": "ns", "qng": "ng",
                   "mudpcu": "pgam", "lambdadpcu": "lamc",
                   "qsde": "qnide"}


def _aero_main():
    m_acc, m_ait, m_crs = AERO_META["mode_idx"]
    l_dst, l_ncl, l_so4 = AERO_META["coarse_species_idx"]
    return dict(scheme="modal",
                nspec=AERO_META["nspec"],
                mode_accum=m_acc - 1, mode_aitken=m_ait - 1,
                mode_coarse=m_crs - 1, coarse_dust=l_dst - 1,
                coarse_nacl=l_ncl - 1, coarse_so4=l_so4 - 1,
                sigmag_aitken=AERO_META["sigmag"][m_ait - 1],
                sigmag_amode=np.array(AERO_META["sigmag"]),
                specdens=D["aero_specdens"],
                spechygro=D["aero_spechygro"],
                voltonumblo=D["aero_voltonumblo"],
                voltonumbhi=D["aero_voltonumbhi"],
                num=D["num_a"], mmr=D["mmr_a"], dgnum=D["dgnum"])


@pytest.fixture(scope="module")
def main_results():
    zc = make_zm_const()
    aero = _aero_main()
    out = {}
    for cfg, first in (("m", True), ("n", False)):
        zp = zm_conv.make_zm_param(
            tau=P["tau"], alfa=P["alfa"], ke=P["ke"],
            dmpdz=P["dmpdz"], tpert_fix=bool(P["tpert_fix"]),
            tpert_fac=P["tpert_fac"], tiedke_add=P["tiedke_add"],
            c0_lnd=P["c0_lnd"], c0_ocn=P["c0_ocn"],
            num_cin=P["num_cin"], limcnv=P["limcnv"],
            mx_bot_lyr_adj=P["mx_bot_lyr_adj"],
            trig_dcape=bool(P["trig_dcape"]),
            trig_ull=bool(P["trig_ull"]),
            clos_dyn_adj=bool(P["clos_dyn_adj"]),
            no_deep_pbl=bool(P["no_deep_pbl"]),
            old_snow=bool(P["old_snow"]),
            zm_microp=True, auto_fac=P["auto_fac"],
            accr_fac=P["accr_fac"], micro_dcs=P["micro_dcs"])
        out[cfg] = zm_conv.zm_conv_main(
            D["t"], D["q"], D["omega"], D["pmid"], D["pint"],
            D["pdel"], D["geos"], D["zmid"], D["zint"], D["pbl_hgt"],
            D["tpert"], D["landfrac"], D["t_star"], D["q_star"],
            float(D["time_step"]), first, zc, zp, aero=aero)
    return out


@pytest.mark.parametrize("cfg", MAIN_CFGS)
def test_main_gather_and_indices(main_results, cfg):
    r = main_results[cfg]
    ng = r["lengath"]
    assert ng == int(D[f"lengath_{cfg}"])
    np.testing.assert_array_equal(
        np.asarray(r["gather_index"][:ng]),
        D[f"gather_index_{cfg}"][:ng] - 1)
    for f in ("msemax_klev_g", "jctop", "jcbot", "jt"):
        np.testing.assert_array_equal(np.asarray(r[f]),
                                      D[f"{f}_{cfg}"] - 1)


# (field, rtol, atol); measured values documented at the bottom of
# the module docstring / PORTING_PLAN row 9
MAIN_REPLAY_SPECS = [
    ("cape", 1e-10, 1e-8),
    ("dcape", 1e-10, 1e-12),
    ("prec", 1e-9, 1e-18),
    ("rliq", 1e-9, 1e-18),
    ("rice", 1e-9, 1e-18),
    ("qtnd", 1e-9, 1e-16),
    ("heat", 2e-9, 1e-12),
    ("mcon", 1e-10, 1e-15),
    ("pflx", 1e-9, 1e-14),
    ("zdu", 1e-10, 1e-16),
    ("mflx_up", 1e-10, 1e-15),
    ("entr_up", 1e-10, 1e-16),
    ("detr_up", 1e-10, 1e-16),
    ("mflx_dn", 1e-10, 1e-15),
    ("entr_dn", 1e-10, 1e-16),
    ("p_del", 1e-13, 0.0),
    ("dsubcld", 1e-13, 0.0),
    ("ql", 1e-9, 1e-13),
    ("rprd", 1e-9, 1e-15),
    ("dlf", 1e-9, 1e-16),
]


@pytest.mark.parametrize("cfg", MAIN_CFGS)
@pytest.mark.parametrize("field,rtol,atol",
                         MAIN_REPLAY_SPECS,
                         ids=[s[0] for s in MAIN_REPLAY_SPECS])
def test_main_replay(main_results, cfg, field, rtol, atol):
    np.testing.assert_allclose(np.asarray(main_results[cfg][field]),
                               D[f"{field}_{cfg}"],
                               rtol=rtol, atol=atol)


@pytest.mark.parametrize("cfg", MAIN_CFGS)
@pytest.mark.parametrize("i", range(len(MICROP_FIELDS)))
def test_main_microp_replay(main_results, cfg, i):
    nm = MICROP_FIELDS[i]
    key = MICROP_NAME_MAP.get(nm, nm)
    got = np.asarray(main_results[cfg]["microp"][key])
    want = D[f"microp_{cfg}"][:, :, i]
    rtol, atol = MTOL.get(nm, MTOL["default"])
    _cmp(f"microp.{nm}", got, want, rtol, atol)


MTOL = dict(KTOL)
MTOL["default"] = (1e-9, 1e-20)
# through the full main (2 outer x 2 inner iterations + closure
# scaling) the bergnn/trspcm noise grows slightly: measured 1.2e-11
# rel (bergnn, values O(10)) and 2.7e-18 abs on trspcm entries of
# O(1e-9) kg/kg/s (the mu*qc(k)-mu*qc(k-1) cancellation again)
MTOL["bergnn"] = (1e-10, 0.0)
MTOL["trspcm"] = (5e-9, 1e-17)


# ---------------------------------------------------------------------------
# Tier-0 properties (microp path)
# ---------------------------------------------------------------------------
class TestTier0:
    def test_warm_family_is_liquid_only(self, kernel_result):
        """Whole-updraft-above-freezing plumes make no ice, snow,
        graupel or freezing heating."""
        r = kernel_result
        w = slice(0, 8)
        for nm in ("qi", "ni", "qni", "ns", "qg", "ng", "sprd", "frz"):
            assert np.all(r[nm][w] == 0.0), nm
        assert r["qc"][w].max() > 0.0 and r["rprd"][w].max() > 0.0

    def test_positivity(self, kernel_result):
        r = kernel_result
        for nm in ("qc", "qi", "nc", "ni", "qr", "nr", "qni", "ns",
                   "qg", "ng", "qcde", "qide", "qnide", "ncde",
                   "nide", "nsde", "rprd", "sprd", "wu"):
            assert np.all(r[nm] >= 0.0), nm

    def test_noop_columns(self, kernel_result):
        """eps0=0 and cmel=cmei=0 columns are exact no-ops."""
        r = kernel_result
        for col in (36, 39):   # dynamics family j=4 (eps0=0), j=7 (cond=0)
            for nm in STATE_FIELDS[:20]:
                assert np.all(r[nm][col] == 0.0), (col, nm)

    def test_aerosol_monotonicity(self, kernel_result):
        """More aerosol -> more activated droplets (aerosol sweep
        family, identical dynamics)."""
        nc = kernel_result["nc"][24:32].max(axis=1)
        assert nc[-1] > nc[0]

    def test_freezing_only_below_freezing(self, kernel_result):
        """Freezing heating never occurs where the whole layer (and
        the one below, where frz is stored) is above 0C."""
        r = kernel_result
        GRAV, CP = ZC["grav"], ZC["cpair"]
        tu = D["k_su"] - GRAV / CP * D["k_zf"][:, :72]
        warm = (tu > 274.0) & (np.roll(tu, -1, axis=1) > 274.0)
        assert np.all(np.abs(r["frz"][warm]) == 0.0)

    def test_number_mass_consistency(self, kernel_result):
        """Zero mass implies zero number for every hydrometeor."""
        r = kernel_result
        for qn, nn in (("qc", "nc"), ("qi", "ni"), ("qr", "nr"),
                       ("qni", "ns"), ("qg", "ng")):
            assert np.all(r[nn][r[qn] == 0.0] == 0.0), (qn, nn)

    def test_main_water_closure(self, main_results):
        """prec + rliq == -column integral of qtnd (the microp path
        conserves water through detrainment of liquid+ice+snow)."""
        for cfg in MAIN_CFGS:
            r = main_results[cfg]
            prec = np.asarray(r["prec"])
            rliq = np.asarray(r["rliq"])
            qint = (np.asarray(r["qtnd"]) * D["pdel"]).sum(axis=1) \
                / ZC["grav"] / 1000.0
            assert np.abs(prec + rliq + qint).max() < 1e-17

    def test_main_snow_le_precip_production(self, main_results):
        """Column-integrated snow production cannot exceed total
        precip production (the pflxs<=pflx fixer guarantees it)."""
        for cfg in MAIN_CFGS:
            r = main_results[cfg]
            sprd = np.asarray(r["microp"]["sprd"])
            rprd = np.asarray(r["rprd"])
            pd = D["pdel"]
            si = (sprd * pd).sum(axis=1)
            ri = (rprd * pd).sum(axis=1)
            assert np.all(si <= ri * (1 + 1e-9) + 1e-12)

    def test_main_rice_matches_detrained_ice(self, main_results):
        for cfg in MAIN_CFGS:
            r = main_results[cfg]
            dif = np.asarray(r["microp"]["dif"])
            dsf = np.asarray(r["microp"]["dsf"])
            rice = ((dif + dsf) * D["pdel"]).sum(axis=1) \
                / ZC["grav"] / 1000.0
            np.testing.assert_allclose(np.asarray(r["rice"]), rice,
                                       rtol=1e-12, atol=1e-20)
