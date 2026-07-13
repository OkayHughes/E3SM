"""Golden replay + Tier-0 property tests for the EAM P3 stratiform
microphysics port (eam_jax/p3/).

Golden: harness/gen_p3_golden.py run in the scream-dev container
against the REAL micro_p3.F90 (eam variant) + micro_p3_utils +
wv_sat_scream, table p3_lookup_table_1.dat-v4.1.1 (see the metadata
note on 4.1.2). Replay drives eam_jax.p3.main.p3_main with the
golden-stored rain tables (test_rain_table_recompute asserts our
regeneration matches them) and the ice table read from the same file.

VALIDATION STRUCTURE (measured, see the bisection notes below):

1. Process rates (p3_tend_out slots 2-35) replay at <= 8.6e-8 rel
   (worst: slot 26 ni_sublim_tend, a ratio of leading-edge tinies;
   all other slots <= 4.1e-10). This validates part1+part2 tightly.
2. Sedimentation sub-kernels replay their own Fortran goldens
   (drv_cloud/rain/ice_sed on identical inputs) at <= 5e-12 rel.
3. The full-chain state comparison is limited by KNIFE-EDGE CFL
   AMPLIFICATION: sedimentation's substep count nstep=int(Co_max+1)
   and the qsmall gates at sedimentation fronts turn the sub-1e-11
   part2 noise (libm tgamma/pow across container-glibc vs XLA) into
   O(1e-3) local differences. Measured proof that this is not a port
   bug: feeding the FORTRAN sedimentation the JAX part2 state
   reproduces the golden mismatch bit-for-bit in magnitude (max rel
   9.05e-4 / 1.34e-3 for a/0 qi/ni_sed — identical to the JAX chain),
   while JAX-vs-Fortran sedimentation on identical inputs agrees to
   4e-13. The chain-amplified fields therefore use a documented
   ENVELOPE: a bounded fraction of points beyond rtol 1e-9 (15-45%
   depending on how much of the ice band a substep flip shifts;
   measured maxima at each assertion) and the worst absolute error
   bounded by a per-field fraction of the field's magnitude
   (measured max ratios are half the bounds used).

Strict fields target 1e-10 rel (measured <= 1.7e-11): the cloud
(liquid) path is deterministic because cloud sedimentation rarely
hits a substep-count tie in these regimes.
"""
import json
from pathlib import Path

import numpy as np
import pytest

from eam_jax.p3 import DEFAULT_OPTS
from eam_jax.p3 import constants as c
from eam_jax.p3 import tables as p3_tables
from eam_jax.p3.main import p3_main
from eam_jax.p3.saturation import qv_sat
from eam_jax.p3.sedimentation import (
    cloud_sedimentation,
    ice_sedimentation,
    rain_sedimentation,
)

GOLD = Path(__file__).resolve().parents[1] / "golden" / "p3_golden.npz"
D = np.load(GOLD)
META = json.loads(str(D["__metadata__"]))
STATE_FIELDS = META["state_fields"]
DIAG_FIELDS = META["diag_fields"]
FLUX_FIELDS = META["flux_fields"]
P = META["params"]
DT = float(META["dt"])
NCOL, NLEV = META["ncol"], META["nlev"]

# host-side table file (same file the Fortran read in the container)
TABLE_DIR = Path(__file__).resolve().parents[3] \
    / "e3sm-inputdata/atm/scream/tables"

OPTS = dict(DEFAULT_OPTS)
assert OPTS == {  # golden params must be the EAMv3 defaults we ship
    "p3_autocon_coeff": P["autocon_coeff"],
    "p3_qc_autocon_expon": P["qc_autocon_expon"],
    "p3_nc_autocon_expon": P["nc_autocon_expon"],
    "p3_accret_coeff": P["accret_coeff"],
    "p3_qc_accret_expon": P["qc_accret_expon"],
    "p3_wbf_coeff": P["wbf_coeff"],
    "p3_mincdnc": P["mincdnc"],
    "p3_max_mean_rain_size": P["max_mean_rain_size"],
    "p3_embryonic_rain_size": P["embryonic_rain_size"],
    "nccnst": P["nccnst"],
}

ICE_TBL, COLL_TBL = p3_tables.read_ice_lookup_tables(
    str(TABLE_DIR / f"p3_lookup_table_1.dat-v{META['table_version']}"),
    META["table_version"])
# replay uses the Fortran-generated rain tables from the golden;
# test_rain_table_recompute asserts our regeneration matches them
TABLES = {
    "ice_table_vals": ICE_TBL,
    "collect_table_vals": COLL_TBL,
    "vn_table_vals": D["vn_table"],
    "vm_table_vals": D["vm_table"],
    "revap_table_vals": D["revap_table"],
    "mu_r_table_vals": D["mu_r_table"],
}

COMMON = dict(
    pres=D["pres"], dz=D["dz"], nc_nuceat_tend=D["nc_nuceat_tend"],
    nccn_prescribed=D["nccn_prescribed"], ni_activated=D["ni_activated"],
    frzimm=D["frzimm"], frzcnt=D["frzcnt"], frzdep=D["frzdep"],
    inv_qc_relvar=D["inv_qc_relvar"], dpres=D["dpres"], exner=D["exner"],
    cld_frac_r=D["cld_frac_r"], cld_frac_l=D["cld_frac_l"],
    cld_frac_i=D["cld_frac_i"])


def run_step(cfg, step):
    cf = META["configs"][cfg]
    st = D[f"{cfg}_in_state_{step}"]
    kw = {k: st[:, :, i] for i, k in enumerate(STATE_FIELDS)}
    out = p3_main(
        DT,
        do_predict_nc=cf["do_predict_nc"],
        do_prescribed_ccn=cf["do_prescribed_ccn"],
        do_precip_off=False,
        use_hetfrz_classnuc=cf["use_hetfrz_classnuc"],
        do_cooper=cf["do_cooper"],
        qv_prev=D[f"{cfg}_in_qvprev_{step}"],
        t_prev=D[f"{cfg}_in_tprev_{step}"],
        tables=TABLES, opts=OPTS, **kw, **COMMON)
    return {k: np.asarray(v) for k, v in out.items()}


_cache = {}


def step_out(cfg, step):
    key = (cfg, step)
    if key not in _cache:
        _cache[key] = run_step(cfg, step)
    return _cache[key]


def check(a, b, rtol, atol=0.0, msg=""):
    a, b = np.asarray(a), np.asarray(b)
    ok = np.isclose(a, b, rtol=rtol, atol=atol)
    if not np.all(ok):
        bad = ~ok
        i = np.unravel_index(np.argmax(np.abs(a - b) * bad), a.shape)
        rel = abs(a[i] - b[i]) / max(abs(b[i]), 1e-300)
        raise AssertionError(
            f"{msg}: {bad.sum()}/{a.size} mismatches; worst at {i}: "
            f"jax={a[i]:.16e} f90={b[i]:.16e} "
            f"abs={abs(a[i]-b[i]):.3e} rel={rel:.3e}")


def check_envelope(a, b, env, msg, frac=0.15, rtol=1e-9, atol=0.0):
    """Knife-edge envelope: >= (1-frac) of the points must match at
    rtol, and the worst absolute error must stay below
    env * max|ref| + atol (atol guards near-zero reference fields).
    The magnitude bound carries the physics; the fraction documents
    the localization of the substep-count flips."""
    a, b = np.asarray(a), np.asarray(b)
    err = np.abs(a - b)
    scale = max(np.abs(b).max(), 1e-30)
    with np.errstate(invalid="ignore"):
        good = err <= rtol * np.abs(b) + max(1e-25 * scale, atol)
    fbad = 1.0 - good.mean()
    assert fbad <= frac, \
        f"{msg}: {100*fbad:.2f}% of points beyond rtol {rtol}"
    assert err.max() <= env * scale + atol, \
        f"{msg}: max abs err {err.max():.3e} > {env:.0e} * scale {scale:.3e}"


ALL_STEPS = [(cfg, s) for cfg, cf in META["configs"].items()
             for s in range(cf["nsteps"])]

# strict fields: measured worst rel <= 1.7e-11 across all configs/steps
# (cancellation atol floors sit at the denormal residue of tendency
# subtractions, measured <= 2e-20 except the interface fluxes)
STRICT_STATE = {
    "qc": (1e-10, 1e-22), "nc": (1e-10, 1e-6), "qv": (1e-10, 0.0),
}
STRICT_DIAG = {
    "diag_eff_radius_qc": (1e-10, 1e-19),
    "qv2qi_depos_tend": (1e-10, 1e-19),
    "precip_total_tend": (1e-10, 1e-19),
    "nevapr": (1e-10, 1e-19),
    "qr_evap_tend": (1e-10, 1e-19),
    "mu_c": (1e-10, 1e-13),
    "lamc": (1e-10, 1e-7),
    "liq_ice_exchange": (1e-10, 1e-19),
    "vap_liq_exchange": (1e-10, 1e-19),
    "vap_ice_exchange": (1e-10, 1e-19),
}
STRICT_FLUX = {
    "precip_liq_flux": (1e-10, 1e-16),
    "precip_ice_flux": (0.0, 0.0),   # identically zero on both sides
    "rflx": (1e-10, 1e-16),
    "cflx": (1e-10, 1e-18),
}

# envelope fields: (max_abs_err / max|ref| bound); measured maxima are
# at most half of these (worst offenders: qr 4.1e-2 c/0, rho_qi 5.6e-2
# a/2, diag_ze_ice 5.5e-2 a/1, sflx 1.6e-2 c/0, th 1.1e-4 b/0)
ENV_STATE = {
    "qr": 8e-2, "nr": 5e-4, "th": 3e-4, "qi": 2e-2, "qm": 2e-2,
    "ni": 2e-3, "bm": 2e-2,
}
ENV_DIAG = {
    "diag_eff_radius_qi": 8e-2, "rho_qi": 1.2e-1,
    "diag_equiv_reflectivity": 6e-2, "diag_ze_rain": 4e-2,
    "diag_ze_ice": 1.2e-1,
}
ENV_FLUX = {"sflx": 4e-2}


@pytest.mark.parametrize("cfg,step", ALL_STEPS)
@pytest.mark.parametrize("field", STATE_FIELDS)
def test_replay_state(cfg, step, field):
    out = step_out(cfg, step)
    ref = D[f"{cfg}_state_{step}"][:, :, STATE_FIELDS.index(field)]
    if field in STRICT_STATE:
        rtol, atol = STRICT_STATE[field]
        check(out[field], ref, rtol, atol, f"{cfg}/{step}/{field}")
    else:
        check_envelope(out[field], ref, ENV_STATE[field],
                       f"{cfg}/{step}/{field}")


@pytest.mark.parametrize("cfg,step", ALL_STEPS)
@pytest.mark.parametrize("field", DIAG_FIELDS)
def test_replay_diag(cfg, step, field):
    out = step_out(cfg, step)
    ref = D[f"{cfg}_diag_{step}"][:, :, DIAG_FIELDS.index(field)]
    if field in STRICT_DIAG:
        rtol, atol = STRICT_DIAG[field]
        check(out[field], ref, rtol, atol, f"{cfg}/{step}/{field}")
    else:
        # measured point fractions reach 19% (sflx b/0) / 16%
        # (reflectivity b/0) where a whole ice band shifts
        check_envelope(out[field], ref, ENV_DIAG[field],
                       f"{cfg}/{step}/{field}", frac=0.30)


@pytest.mark.parametrize("cfg,step", ALL_STEPS)
@pytest.mark.parametrize("field", FLUX_FIELDS)
def test_replay_flux(cfg, step, field):
    out = step_out(cfg, step)
    ref = D[f"{cfg}_flux_{step}"][:, :, FLUX_FIELDS.index(field)]
    if field in STRICT_FLUX:
        rtol, atol = STRICT_FLUX[field]
        check(out[field], ref, rtol, atol, f"{cfg}/{step}/{field}")
    else:
        check_envelope(out[field], ref, ENV_FLUX[field],
                       f"{cfg}/{step}/{field}", frac=0.30)


@pytest.mark.parametrize("cfg,step", ALL_STEPS)
def test_replay_surf(cfg, step):
    out = step_out(cfg, step)
    ref = D[f"{cfg}_surf_{step}"]
    check(out["precip_liq_surf"], ref[:, 0], 1e-10, 1e-19,
          f"{cfg}/{step}/precip_liq_surf")
    # ice surface precip inherits the ice-sed knife edge in EVERY
    # ice-producing column (measured: 35% of the 40 columns beyond
    # 1e-9 at b/0, max abs/scale 3.5e-2)
    check_envelope(out["precip_ice_surf"], ref[:, 1], 8e-2,
                   f"{cfg}/{step}/precip_ice_surf", frac=0.45)


@pytest.mark.parametrize("cfg,step", ALL_STEPS)
def test_replay_process_rates(cfg, step):
    """p3_tend_out slots 2-35 (the part2 process rates) replay tightly:
    measured worst 8.6e-8 rel (slot 26 ni_sublim_tend at a/2 — the
    ni_incld/qi_incld ratio of leading-edge tinies; all other slots
    <= 4.1e-10)."""
    out = step_out(cfg, step)
    ref = D[f"{cfg}_tend_{step}"]
    got = out["p3_tend_out"]
    check(got[..., 1:35], ref[..., 1:35], 5e-7, 1e-19,
          f"{cfg}/{step}/tend[2-35]")


@pytest.mark.parametrize("cfg,step", ALL_STEPS)
def test_replay_sed_and_measured_tends(cfg, step):
    """Slots 36-41 (sedimentation) and 42-49 (measured state deltas)
    carry the knife-edge envelope and, for the deltas, the (a-b)/dt
    cancellation noise (rel 1.0 on true-zero tendencies)."""
    out = step_out(cfg, step)
    ref = D[f"{cfg}_tend_{step}"]
    got = out["p3_tend_out"]
    for slot in range(36, 50):
        # atol floor: slot 37 (nc sed) under predict_nc=F is pure
        # dsd-consistency noise with a ~1e-11 #/kg/s scale
        check_envelope(got[..., slot - 1], ref[..., slot - 1], 5e-2,
                       f"{cfg}/{step}/tend[{slot}]", frac=0.30,
                       atol=1e-10)


# ------------------------------------------------------------------
# sedimentation sub-kernel replay (tight; inputs recorded in golden)
# ------------------------------------------------------------------

def _sed_common():
    return dict(rho=D["sed_rho"], inv_rho=D["sed_inv_rho"],
                inv_dz=1.0 / D["dz"])


def test_sed_kernel_cloud():
    s = _sed_common()
    st0 = D["a_in_state_0"]
    qc = st0[:, :, STATE_FIELDS.index("qc")]
    nc = st0[:, :, STATE_FIELDS.index("nc")]
    zflx = np.zeros((NCOL, NLEV + 1))
    zc = np.zeros((NCOL, NLEV))
    r = cloud_sedimentation(
        D["sed_qc_incld"], s["rho"], s["inv_rho"], D["cld_frac_l"],
        D["sed_acn"], s["inv_dz"], DT, 1.0 / DT, True,
        qc, nc, D["sed_nc_incld"], zc, zc,
        np.zeros(NCOL), zflx, qc, nc)
    for k, gk in [("qc", "csed_qc"), ("nc", "csed_nc"),
                  ("nc_incld", "csed_nc_incld"), ("mu_c", "csed_mu_c"),
                  ("lamc", "csed_lamc"), ("cflx", "csed_cflx"),
                  ("qc_tend", "csed_qc_tend"), ("nc_tend", "csed_nc_tend"),
                  ("precip_liq_surf", "csed_precip_liq_surf")]:
        check(np.asarray(r[k]), D[gk], 5e-12, 1e-20, f"cloud sed {k}")


def test_sed_kernel_rain():
    s = _sed_common()
    st0 = D["a_in_state_0"]
    qr = st0[:, :, STATE_FIELDS.index("qr")]
    nr = st0[:, :, STATE_FIELDS.index("nr")]
    zflx = np.zeros((NCOL, NLEV + 1))
    zc = np.zeros((NCOL, NLEV))
    r = rain_sedimentation(
        s["rho"], s["inv_rho"], D["sed_rhofacr"], D["cld_frac_r"],
        s["inv_dz"], D["sed_qr_incld"],
        TABLES["vn_table_vals"], TABLES["vm_table_vals"], DT, 1.0 / DT,
        OPTS["p3_max_mean_rain_size"], qr, nr, D["sed_nr_incld"], zc, zc,
        np.zeros(NCOL), zflx, zflx, qr, nr)
    for k, gk in [("qr", "rsed_qr"), ("nr", "rsed_nr"),
                  ("nr_incld", "rsed_nr_incld"), ("mu_r", "rsed_mu_r"),
                  ("lamr", "rsed_lamr"),
                  ("precip_liq_flux", "rsed_precip_liq_flux"),
                  ("rflx", "rsed_rflx"), ("qr_tend", "rsed_qr_tend"),
                  ("nr_tend", "rsed_nr_tend"),
                  ("precip_liq_surf", "rsed_precip_liq_surf")]:
        check(np.asarray(r[k]), D[gk], 5e-12, 1e-20, f"rain sed {k}")


def test_sed_kernel_ice():
    s = _sed_common()
    st0 = D["a_in_state_0"]
    qi = st0[:, :, STATE_FIELDS.index("qi")]
    ni = st0[:, :, STATE_FIELDS.index("ni")]
    qm = st0[:, :, STATE_FIELDS.index("qm")]
    bm = st0[:, :, STATE_FIELDS.index("bm")]
    zflx = np.zeros((NCOL, NLEV + 1))
    r = ice_sedimentation(
        s["rho"], s["inv_rho"], D["sed_rhofaci"], D["cld_frac_i"],
        s["inv_dz"], DT, 1.0 / DT,
        qi, D["sed_qi_incld"], ni, D["sed_ni_incld"],
        qm, D["sed_qm_incld"], bm, D["sed_bm_incld"],
        TABLES["ice_table_vals"], np.zeros(NCOL), zflx, qi, ni)
    for k, gk in [("qi", "ised_qi"), ("ni", "ised_ni"), ("qm", "ised_qm"),
                  ("bm", "ised_bm"), ("qi_incld", "ised_qi_incld"),
                  ("ni_incld", "ised_ni_incld"),
                  ("qm_incld", "ised_qm_incld"),
                  ("bm_incld", "ised_bm_incld"), ("sflx", "ised_sflx"),
                  ("qi_tend", "ised_qi_tend"), ("ni_tend", "ised_ni_tend"),
                  ("precip_ice_surf", "ised_precip_ice_surf")]:
        check(np.asarray(r[k]), D[gk], 5e-12, 1e-18, f"ice sed {k}")


def test_qv_sat_probe():
    """Unit bisection of the saturation port (both phases)."""
    got_l = np.asarray(qv_sat(D["qsat_probe_t"], D["qsat_probe_p"], False))
    got_i = np.asarray(qv_sat(D["qsat_probe_t"], D["qsat_probe_p"], True))
    check(got_l, D["qsat_probe_liq"], 1e-13, 0.0, "qv_sat liq")
    check(got_i, D["qsat_probe_ice"], 1e-13, 0.0, "qv_sat ice")


def test_rain_table_recompute():
    """p3_init_b regeneration matches the Fortran tables (summation
    order in the 10000-bin PSD integral is the only difference)."""
    mu_r, vn, vm, revap = p3_tables.compute_rain_tables()
    check(mu_r, D["mu_r_table"], 0.0, 0.0, "mu_r table")
    check(vn, D["vn_table"], 1e-12, 0.0, "vn table")
    check(vm, D["vm_table"], 1e-12, 0.0, "vm table")
    check(revap, D["revap_table"], 1e-12, 0.0, "revap table")


# ------------------------------------------------------------------
# Tier-0 property tests
# ------------------------------------------------------------------

def _column_water(state_arr, dpres):
    idx = [STATE_FIELDS.index(k) for k in
           ("qv", "qc", "qr", "qi")]  # qm is part of qi, not extra mass
    tot = sum(state_arr[:, :, i] for i in idx)
    return np.sum(tot * dpres / c.g, axis=-1)


@pytest.mark.parametrize("cfg,step", ALL_STEPS)
def test_tier0_water_conservation(cfg, step):
    """Column water in == column water out + surface precip mass."""
    out = step_out(cfg, step)
    st_in = D[f"{cfg}_in_state_{step}"].copy()
    # p3_main clips negative qv on entry, exactly like the Fortran
    st_in[:, :, STATE_FIELDS.index("qv")] = \
        np.maximum(st_in[:, :, STATE_FIELDS.index("qv")], 0.0)
    w_in = _column_water(st_in, D["dpres"])
    st_out = np.stack([out[k] for k in STATE_FIELDS], axis=-1)
    w_out = _column_water(st_out, D["dpres"])
    precip = (out["precip_liq_surf"] + out["precip_ice_surf"]) \
        * c.rho_h2o * DT
    resid = w_in - (w_out + precip)
    assert np.max(np.abs(resid) / np.maximum(w_in, 1.0)) < 1e-12, \
        f"water closure residual {np.max(np.abs(resid)):.3e} kg/m2"


def test_tier0_positivity():
    for cfg, step in ALL_STEPS:
        out = step_out(cfg, step)
        for f in ("qc", "qr", "qi", "qm", "qv", "bm"):
            assert out[f].min() >= 0.0, f"{cfg}/{step}: negative {f}"
        assert out["precip_liq_surf"].min() >= 0.0
        assert out["precip_ice_surf"].min() >= 0.0


def test_tier0_cloud_free_noop():
    """The fully clear column (family 4, j=0 -> column 32) must pass
    through untouched: no tendencies, no precip, default diagnostics."""
    i = 32
    out = step_out("a", 0)
    st_in = D["a_in_state_0"]
    for k in ("qc", "nc", "qr", "nr", "qi", "ni", "qm", "bm"):
        j = STATE_FIELDS.index(k)
        np.testing.assert_array_equal(out[k][i], st_in[i, :, j])
    assert out["precip_liq_surf"][i] == 0.0
    assert out["precip_ice_surf"][i] == 0.0
    assert np.all(out["p3_tend_out"][i] == 0.0)
    np.testing.assert_array_equal(out["diag_eff_radius_qc"][i], 10.0e-6)
    np.testing.assert_array_equal(out["diag_equiv_reflectivity"][i], -99.0)


def test_tier0_nucleation_only_column():
    """Cold clear ice-supersaturated column (33) must gain ice from
    deposition nucleation alone."""
    out = step_out("a", 0)
    assert out["qi"][33].max() > 0.0
    assert out["ni"][33].max() > 0.0


def test_tier0_energy_consistency():
    """Column liquid-water-ice static energy analogue
    cp*T - Lv*(qc+qr) - Ls*qi closes against the precip fluxes."""
    out = step_out("a", 0)
    st_in = D["a_in_state_0"].copy()
    st_in[:, :, STATE_FIELDS.index("qv")] = \
        np.maximum(st_in[:, :, STATE_FIELDS.index("qv")], 0.0)
    exner = D["exner"]
    dp_g = D["dpres"] / c.g

    def energy(th, qc, qr, qi):
        t = th / exner
        return np.sum((c.cp * t - c.latvap * (qc + qr)
                       - c.latsub * qi) * dp_g, axis=-1)

    idx = {k: STATE_FIELDS.index(k) for k in STATE_FIELDS}
    e_in = energy(st_in[:, :, idx["th"]], st_in[:, :, idx["qc"]],
                  st_in[:, :, idx["qr"]], st_in[:, :, idx["qi"]])
    e_out = energy(out["th"], out["qc"], out["qr"], out["qi"])
    # precip leaves with its latent-energy deficit
    e_out = e_out - c.latvap * out["precip_liq_surf"] * c.rho_h2o * DT \
        - c.latsub * out["precip_ice_surf"] * c.rho_h2o * DT
    resid = np.abs(e_in - e_out) / np.maximum(np.abs(e_in), 1.0)
    assert resid.max() < 1e-10, f"energy residual {resid.max():.3e}"
