"""Tier-1: scream_jax RRTMGP vs the EAMxx golden archive.

Replays golden/rrtmgp_218x72_dt1800_3steps.npz: starting from each
captured state n-1, runs one rrtmgp_process_step (radiation updates
every step in the archive) and compares every computed EAMxx field
against the captured state n. Column lat/lon comes from the IC file
recorded in the metadata.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

REPO = Path(__file__).resolve().parents[2]
GOLDEN = REPO / "jax_port" / "golden" / "rrtmgp_218x72_dt1800_3steps.npz"
DDIR = REPO.parent / "e3sm-inputdata" / "atm" / "scream" / "init"

from scream_jax.rrtmgp.coefficients import load_kdist  # noqa: E402
from scream_jax.rrtmgp.cloud_optics import load_cloud_optics  # noqa: E402
from scream_jax.rrtmgp.process import (DEFAULT_PARAMS, GAS_NAMES,  # noqa: E402
                                       rrtmgp_process_step)

# state fed forward between steps
PROGNOSTICS = ("T_mid", "rad_heating_pdel")
FORCINGS = ("p_mid", "p_int", "pseudo_density", "qv", "qc", "nc", "qi",
            "cldfrac_tot", "eff_radius_qc", "eff_radius_qi",
            "sfc_alb_dir_vis", "sfc_alb_dir_nir", "sfc_alb_dif_vis",
            "sfc_alb_dif_nir", "surf_lw_flux_up", "o3_volume_mix_ratio")

REL_TOL = 1e-6
MAX_BAD_FRACTION = 1e-3


def _load():
    if not GOLDEN.exists():
        pytest.skip(f"golden archive not found: {GOLDEN}")
    if not (DDIR / "rrtmgp-data-sw-g112-210809.nc").exists():
        pytest.skip("rrtmgp coefficient files not available")
    z = np.load(GOLDEN)
    meta = json.loads(str(z["__metadata__"]))
    ic = meta["ic_file"]
    cands = [Path(ic)]
    if ic.startswith("/work/e3sm-inputdata"):
        cands.append(REPO.parent / "e3sm-inputdata"
                     / Path(ic).relative_to("/work/e3sm-inputdata"))
    ic_path = next((p for p in cands if p.exists()), None)
    if ic_path is None:
        pytest.skip(f"IC file not found: {ic}")
    nc4 = pytest.importorskip("netCDF4")
    ds = nc4.Dataset(ic_path)
    lat = np.array(ds["lat"][:])
    lon = np.array(ds["lon"][:])
    ds.close()
    return z, meta, lat, lon


def test_rrtmgp_golden_replay():
    z, meta, lat, lon = _load()
    dt = float(meta["dt"])
    # t0 2021-10-12-45000: day-of-year 285, plus seconds fraction; the
    # C++ uses calday = frac_of_year_in_days() + 1
    doy0 = 285
    sec0 = 45000.0
    year = 2021

    kd_sw = load_kdist(str(DDIR / "rrtmgp-data-sw-g112-210809.nc"), GAS_NAMES)
    kd_lw = load_kdist(str(DDIR / "rrtmgp-data-lw-g128-210809.nc"), GAS_NAMES)
    co_sw = load_cloud_optics(str(DDIR / "rrtmgp-cloud-optics-coeffs-sw.nc"))
    co_lw = load_cloud_optics(str(DDIR / "rrtmgp-cloud-optics-coeffs-lw.nc"))

    params = dict(DEFAULT_PARAMS)
    params["orbital_year"] = 1990
    params["rad_frequency"] = 1

    def f(name, step):
        return z[f"{name}__step{step}"]

    worst, worst_frac = {}, {}
    for step in range(1, meta["nsteps"] + 1):
        prev = step - 1
        nstep = prev  # start-of-step step count
        calday = (doy0 - 1) + (sec0 + prev * dt) / 86400.0 + 1

        out = rrtmgp_process_step(
            kd_sw, kd_lw, co_sw, co_lw, params,
            dt, nstep, year, calday, lat, lon,
            f("T_mid", prev), f("p_mid", prev), f("p_int", prev),
            f("pseudo_density", prev),
            f("sfc_alb_dir_vis", prev), f("sfc_alb_dir_nir", prev),
            f("sfc_alb_dif_vis", prev), f("sfc_alb_dif_nir", prev),
            f("qv", prev), f("qc", prev), f("nc", prev), f("qi", prev),
            f("cldfrac_tot", prev),
            f("eff_radius_qc", prev), f("eff_radius_qi", prev),
            f("surf_lw_flux_up", prev), f("o3_volume_mix_ratio", prev),
            f("rad_heating_pdel", prev))

        for name in meta["computed_fields"]:
            gname = f"{name}__step{step}"
            if gname not in z.files or name not in out:
                continue
            mine = np.asarray(out[name], dtype=np.float64)
            ref = np.asarray(z[gname], dtype=np.float64)
            scale = max(np.abs(ref).max(), 1e-30)
            rel = np.abs(mine - ref) / scale
            bad = rel > REL_TOL
            worst_frac[name] = max(worst_frac.get(name, 0.0), bad.mean())
            worst[name] = max(worst.get(name, 0.0), rel.max())

    report = "\n".join(
        f"  {k:28s} max rel err = {worst[k]:.3e}  "
        f"(outlier fraction {worst_frac[k]:.2e})"
        for k in sorted(worst))
    print(f"\nTier-1 RRTMGP golden replay ({meta['nsteps']} steps):\n{report}")

    over = {k: v for k, v in worst_frac.items() if v > MAX_BAD_FRACTION}
    assert not over, (
        f"fields with more than {MAX_BAD_FRACTION:.1%} of points beyond "
        f"{REL_TOL}:\n{report}")
