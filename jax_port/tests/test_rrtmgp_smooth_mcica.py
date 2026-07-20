"""Smoothed MCICA subcolumn masks (approximation-by-identity).

The MCICA subcolumn mask in scream_jax.rrtmgp.mcica is a composition of
hard comparisons on the JSF64 random deviates, so d(mask)/d(cldfrac) is
identically zero under AD. The static kwarg smooth_width (threaded
rrtmgp_process_step -> rrtmgp_main -> get_subsampled_clouds ->
get_subcolumn_mask) relaxes those comparisons: width 0 is the exact
bitwise scheme, width > 0 yields fractional masks that multiply the
cloud optics, making d(fluxes)/d(cldfrac_tot) respond through the
sampling.

Tests here check:
- width 0 output is array-equal to a call without the kwarg (bitwise
  default preserved);
- width > 0 fluxes are finite, physically bounded, and converge to the
  hard fluxes as width -> 0;
- the mask path gradient contrast: identically zero at width 0,
  nonzero/finite at width 0.1 (asserted at the get_subsampled_clouds
  level with FIXED optics so only the mask depends on cldfrac), and at
  the process level the width-0.1 gradient is finite, NaN-free, nonzero
  and much larger than the width-0 gradient (which is NOT identically
  zero: a smooth in-cloud-water path through lwp/iwp exists even in
  hard mode);
- grad w.r.t. T_mid at small width stays close to the hard-mode value;
- a central finite-difference check of the SMOOTHED primal against the
  smoothed gradient in a cldfrac direction supported on fractional-
  cloud cells (chosen away from the hard branches that remain outside
  MCICA scope: the cf > 0 jump and the in-cloud min() kink in
  mixing_ratio_to_cloud_mass).

Uses the golden-state loading pattern of test_ad_rrtmgp.py, but from
STEP 1 and from columns that actually contain fractional cloud (step 0
of the golden archive is cloud-free, which makes every cldfrac
derivative trivially zero).
"""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import jax  # noqa: E402
import jax.numpy as jnp  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
GOLDEN = REPO / "jax_port" / "golden" / "physics_suite_218x72_dt1800_2steps.npz"
DATA = REPO.parent / "e3sm-inputdata" / "atm" / "scream"
IC_FILE = DATA / "init" / "screami_unit_tests_ne2np4L72_20220822.nc"

# step-1 columns of the golden archive with fractional cloud cover
COLS = np.array([2, 3, 5, 6, 7, 13])

# widths for the convergence scan, descending
WIDTHS = (0.3, 0.1, 0.03, 0.01, 0.003)

FLUX_KEYS = ("SW_flux_up", "SW_flux_dn", "SW_flux_dn_dir",
             "LW_flux_up", "LW_flux_dn")


@pytest.fixture(scope="module")
def setup():
    """Golden state, physics tables and a rad(cldfrac, T, width) closure
    plus the shared (expensive) primal/gradient evaluations."""
    if not GOLDEN.exists() or not IC_FILE.exists():
        pytest.skip("golden archive or input data not available")
    nc4 = pytest.importorskip("netCDF4")

    z = np.load(GOLDEN)
    meta = json.loads(str(z["__metadata__"]))
    s = {name: np.asarray(z[f"{name}__step1"], dtype=np.float64)[COLS]
         for name in meta["fields"]}
    ds = nc4.Dataset(IC_FILE)
    geo = {k: np.array(ds[k][:]) for k in ("hyam", "hybm", "lat", "lon",
                                           "area")}
    ds.close()

    from scream_jax.driver import ScreamPhysics
    from scream_jax.foundation.thermo import calculate_dx_from_area
    from scream_jax.rrtmgp.process import rrtmgp_process_step

    phys = ScreamPhysics(DATA, geo["hyam"], geo["hybm"], geo["lat"][COLS],
                         geo["lon"][COLS],
                         np.asarray(calculate_dx_from_area(
                             geo["area"][COLS], geo["lat"][COLS])))

    def rad(cldfrac, T, width=None):
        kw = {} if width is None else {"smooth_width": width}
        return rrtmgp_process_step(
            phys.kd_sw, phys.kd_lw, phys.co_sw, phys.co_lw,
            phys.rrtmgp_params, 1800.0, 0, 2021, 285.5,
            jnp.asarray(phys.lat), jnp.asarray(phys.lon),
            T, jnp.asarray(s["p_mid"]), jnp.asarray(s["p_int"]),
            jnp.asarray(s["pseudo_density"]),
            jnp.asarray(s["sfc_alb_dir_vis"]),
            jnp.asarray(s["sfc_alb_dir_nir"]),
            jnp.asarray(s["sfc_alb_dif_vis"]),
            jnp.asarray(s["sfc_alb_dif_nir"]),
            jnp.asarray(s["qv"]), jnp.asarray(s["qc"]),
            jnp.asarray(s["nc"]), jnp.asarray(s["qi"]),
            cldfrac,
            jnp.asarray(s["eff_radius_qc"]),
            jnp.asarray(s["eff_radius_qi"]),
            jnp.asarray(s["surf_lw_flux_up"]),
            jnp.asarray(s["o3_volume_mix_ratio"]),
            jnp.asarray(s["rad_heating_pdel"]),
            aero_tau_sw=jnp.asarray(s["aero_tau_sw"]),
            aero_ssa_sw=jnp.asarray(s["aero_ssa_sw"]),
            aero_g_sw=jnp.asarray(s["aero_g_sw"]),
            aero_tau_lw=jnp.asarray(s["aero_tau_lw"]),
            **kw)

    c0 = jnp.asarray(s["cldfrac_tot"])
    T0 = jnp.asarray(s["T_mid"])

    def scalar_c(c, width):
        out = rad(c, T0, width)
        return jnp.sum(out["SW_flux_dn"] + out["LW_flux_up"])

    def scalar_T(T, width):
        return jnp.sum(rad(c0, T, width)["LW_flux_up"])

    out_default = rad(c0, T0)             # no kwarg
    out_hard = rad(c0, T0, 0.0)           # explicit width 0
    out_w = {w: rad(c0, T0, w) for w in WIDTHS}

    grad_c = {w: np.asarray(jax.grad(lambda c: scalar_c(c, w))(c0))
              for w in (0.0, 0.1)}
    grad_T_hard = np.asarray(jax.grad(lambda T: scalar_T(T, 0.0))(T0))
    grad_T_small = np.asarray(jax.grad(lambda T: scalar_T(T, 0.01))(T0))

    return {"s": s, "phys": phys, "rad": rad, "c0": c0, "T0": T0,
            "scalar_c": scalar_c,
            "out_default": out_default, "out_hard": out_hard,
            "out_w": out_w, "grad_c": grad_c,
            "grad_T_hard": grad_T_hard, "grad_T_small": grad_T_small}


def test_state_has_fractional_cloud(setup):
    """Guard: the chosen columns exercise fractional cloud cover."""
    c = np.asarray(setup["c0"])
    assert int(((c > 0) & (c < 1)).sum()) >= 10


def test_width_zero_matches_no_kwarg(setup):
    """smooth_width=0.0 is bitwise the default (no-kwarg) behavior."""
    for key in FLUX_KEYS + ("T_mid", "rad_heating_pdel", "cldtot"):
        a = np.asarray(setup["out_default"][key])
        b = np.asarray(setup["out_hard"][key])
        assert np.array_equal(a, b), f"{key} differs at width 0"


def test_width_zero_mask_bitwise(setup):
    """get_subcolumn_mask at width 0 equals the no-kwarg mask exactly
    (dtype included)."""
    from scream_jax.rrtmgp import mcica
    c = setup["c0"]
    seeds = mcica.compute_seeds(jnp.asarray(setup["s"]["p_mid"]), 1)
    ngpt = int(setup["phys"].kd_sw["ngpt"])
    m_ref = np.asarray(mcica.get_subcolumn_mask(c, seeds, ngpt))
    m_0 = np.asarray(mcica.get_subcolumn_mask(c, seeds, ngpt,
                                              smooth_width=0.0))
    assert m_ref.dtype == m_0.dtype == np.int32
    assert np.array_equal(m_ref, m_0)


def test_smoothed_fluxes_finite_and_physical(setup):
    """Width > 0 fluxes are finite and within loose physical bounds."""
    bounds = {"SW_flux_up": (0.0, 1500.0), "SW_flux_dn": (0.0, 1500.0),
              "SW_flux_dn_dir": (0.0, 1500.0),
              "LW_flux_up": (0.0, 1000.0), "LW_flux_dn": (0.0, 1000.0)}
    for w, out in setup["out_w"].items():
        for key, (lo, hi) in bounds.items():
            f = np.asarray(out[key])
            assert np.isfinite(f).all(), f"{key} not finite at width {w}"
            assert f.min() >= lo and f.max() <= hi, \
                f"{key} out of [{lo},{hi}] at width {w}: " \
                f"[{f.min()}, {f.max()}]"


def test_smoothed_fluxes_converge_to_hard(setup):
    """max|flux(width) - flux(0)| decreases monotonically as width -> 0.

    Measured trend (ncol=6 golden cols, max over the 5 flux fields):
    width 0.3 -> ~1.0e3 W/m2, 0.1 -> ~7.7e2, 0.03 -> ~4.2e2,
    0.01 -> ~1.6e2, 0.003 -> ~1.3e2."""
    hard = {k: np.asarray(setup["out_hard"][k]) for k in FLUX_KEYS}
    dmax = []
    for w in WIDTHS:  # descending widths
        out = setup["out_w"][w]
        dmax.append(max(float(np.abs(np.asarray(out[k]) - hard[k]).max())
                        for k in FLUX_KEYS))
    assert all(b < a for a, b in zip(dmax, dmax[1:])), \
        f"not monotone: {list(zip(WIDTHS, dmax))}"
    assert dmax[-1] < 0.25 * dmax[0], \
        f"weak convergence: {list(zip(WIDTHS, dmax))}"


def test_mask_path_gradient_contrast_unit_level(setup):
    """With FIXED cloud optics (so cldfrac only enters through the
    mask), d(sum subsampled tau)/d(cldfrac) is identically zero at
    width 0 and nonzero, NaN-free at width 0.1."""
    from scream_jax.rrtmgp import mcica
    phys = setup["phys"]
    c0 = setup["c0"]
    ncol, nlay = np.asarray(c0).shape
    nbnd, ngpt = int(phys.kd_sw["nband"]), int(phys.kd_sw["ngpt"])
    rng = np.random.default_rng(1)
    opt = {"tau": jnp.asarray(rng.uniform(0.1, 5.0, (ncol, nlay, nbnd))),
           "ssa": jnp.asarray(rng.uniform(0.3, 0.9, (ncol, nlay, nbnd))),
           "g": jnp.asarray(rng.uniform(0.1, 0.8, (ncol, nlay, nbnd)))}
    play = jnp.asarray(setup["s"]["p_mid"])

    def tau_sum(c, w):
        o = mcica.get_subsampled_clouds(opt, c, play,
                                        phys.kd_sw["gpt2band"], ngpt,
                                        two_stream=True, smooth_width=w)
        return jnp.sum(o["tau"])

    g0 = np.asarray(jax.grad(lambda c: tau_sum(c, 0.0))(c0))
    g1 = np.asarray(jax.grad(lambda c: tau_sum(c, 0.1))(c0))
    assert np.all(g0 == 0.0), "hard mask path should have zero gradient"
    assert np.isfinite(g1).all()
    assert np.linalg.norm(g1) > 0.0
    # every cell sees the smoothed threshold
    assert int((g1 != 0).sum()) == g1.size


def test_process_gradient_cldfrac_smoothed(setup):
    """Reverse grad of sum(SW_flux_dn + LW_flux_up) w.r.t. cldfrac_tot
    at width 0.1: finite, NaN-free, nonzero, and dominated by the mask
    path (norm far above the width-0 gradient, which is NOT identically
    zero because of the smooth lwp/iwp in-cloud-water path)."""
    g0 = setup["grad_c"][0.0]
    g1 = setup["grad_c"][0.1]
    assert np.isnan(g1).sum() == 0 and np.isinf(g1).sum() == 0
    assert np.linalg.norm(g1) > 0.0
    # contrast with hard mode: the mask channel adds orders of magnitude
    n0, n1 = np.linalg.norm(g0), np.linalg.norm(g1)
    assert n1 > 10.0 * n0, f"|g(0.1)|={n1:.3e} vs |g(0)|={n0:.3e}"
    # measured: |g(0)| ~ 1.5e3 (water path only), |g(0.1)| ~ 1.8e6


def test_grad_T_small_width_matches_hard(setup):
    """d(sum LW_flux_up)/d(T_mid) at width 0.01 stays within a loose
    tolerance of the hard-mode gradient (measured rel diff ~1.6e-2)."""
    gh = setup["grad_T_hard"]
    gs = setup["grad_T_small"]
    assert np.isfinite(gh).all() and np.isfinite(gs).all()
    nh = np.linalg.norm(gh)
    assert nh > 0.0
    rel = np.linalg.norm(gs - gh) / nh
    assert rel < 0.1, f"grad_T rel diff {rel:.3e}"


def test_fd_matches_smoothed_gradient(setup):
    """Central FD of the SMOOTHED primal (width 0.1) vs the smoothed
    reverse gradient, in a cldfrac direction supported on fractional-
    cloud cells away from the remaining (out-of-scope) hard branches:
    the cf > 0 jump and the min(q/cf, 0.005) kink in
    mixing_ratio_to_cloud_mass. The smoothed MCICA path is C-infinity,
    so agreement should be tight (measured rel ~8e-8 at eps 1e-3)."""
    c0 = setup["c0"]
    cnp = np.asarray(c0)
    s = setup["s"]
    safe = (cnp > 0.05) & (cnp < 0.95)
    safe &= np.abs(cnp - s["qc"] / 0.005) > 5e-3
    safe &= np.abs(cnp - s["qi"] / 0.005) > 5e-3
    assert safe.sum() >= 3, "not enough branch-safe fractional cells"
    rng = np.random.default_rng(0)
    v = rng.normal(size=cnp.shape) * safe
    v /= np.linalg.norm(v)

    g1 = setup["grad_c"][0.1]
    tang = float(g1.ravel() @ v.ravel())
    assert tang != 0.0

    eps = 1e-3
    scalar_c = setup["scalar_c"]
    fp = float(scalar_c(c0 + eps * jnp.asarray(v), 0.1))
    fm = float(scalar_c(c0 - eps * jnp.asarray(v), 0.1))
    fd = (fp - fm) / (2.0 * eps)
    rel = abs(fd - tang) / max(abs(tang), 1e-30)
    assert rel < 1e-3, f"fd={fd!r} grad.v={tang!r} rel={rel:.3e}"
