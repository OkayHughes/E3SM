"""Tier-1: scream_jax SHOC vs the EAMxx golden archive.

Replays jax_port/golden/shoc_218x72_dt1800_5steps.npz: starting from each
captured state n-1, runs 6 AD subcycles (dt=300 s) of shoc_process_step and
compares every mapped EAMxx field against the captured state n.

Grid geometry (hyam/hybm for npbl, area/lat for cell_length) comes from the
IC file recorded in the archive metadata; the test is skipped if either
file is unavailable.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

REPO = Path(__file__).resolve().parents[2]
GOLDEN = REPO / "jax_port" / "golden" / "shoc_218x72_dt1800_5steps.npz"

from scream_jax.foundation import constants as c  # noqa: E402
from scream_jax.foundation.thermo import calculate_dx_from_area  # noqa: E402
from scream_jax.shoc import shoc_init  # noqa: E402
from scream_jax.shoc.process import shoc_process_step  # noqa: E402

# Runtime params used by gen_golden.py (must match PROC_CONFIGS["shoc"])
PARAMS = dict(lambda_low=0.001, lambda_high=0.08, lambda_slope=2.65,
              lambda_thresh=0.02, thl2tune=1.0, qw2tune=1.0, qwthl2tune=1.0,
              w2tune=1.0, length_fac=0.5, c_diag_3rd_mom=7.0, ckh=0.1, ckm=0.1)
N_SUBCYCLES = 6

# EAMxx field name -> shoc_process_step output key (state fields advanced in time)
PROGNOSTIC_MAP = {
    "T_mid": "T_mid", "qv": "qv", "qc": "qc", "tke": "tke",
    "cldfrac_liq": "cldfrac_liq", "sgs_buoy_flux": "sgs_buoy_flux",
    "eddy_diff_mom": "eddy_diff_mom",
}
DIAGNOSTIC_MAP = {
    "eddy_diff_heat": "eddy_diff_heat", "inv_qc_relvar": "inv_qc_relvar",
    "pbl_height": "pbl_height", "ustar": "ustar", "obklen": "obklen",
    "w_variance": "w_variance", "thl_sec": "thl_sec",
    "cldfrac_liq_prev": "cldfrac_liq_prev",
}


def _load():
    if not GOLDEN.exists():
        pytest.skip(f"golden archive not found: {GOLDEN}")
    z = np.load(GOLDEN)
    meta = json.loads(str(z["__metadata__"]))
    ic_container = meta["ic_file"]
    candidates = [
        Path(ic_container),
        REPO.parent / "e3sm-inputdata" / Path(ic_container).relative_to("/work/e3sm-inputdata")
        if ic_container.startswith("/work/e3sm-inputdata") else None,
    ]
    ic_path = next((p for p in candidates if p and p.exists()), None)
    if ic_path is None:
        pytest.skip(f"IC file not found on host: {ic_container}")
    nc = pytest.importorskip("netCDF4")
    ds = nc.Dataset(ic_path)
    geo = {k: np.array(ds[k][:]) for k in ("hyam", "hybm", "area", "lat")}
    ds.close()
    return z, meta, geo


def test_shoc_golden_replay():
    z, meta, geo = _load()
    dt_host = meta["dt"]
    dt = dt_host / N_SUBCYCLES
    assert meta["params"]["number_of_subcycles"] == N_SUBCYCLES

    pref_mid = c.P0 * (geo["hyam"] + geo["hybm"])
    npbl = shoc_init(meta["nlevs"], 0, pref_mid) if "nlevs" in meta else \
        shoc_init(72, 0, pref_mid)
    cell_length = np.asarray(calculate_dx_from_area(geo["area"], geo["lat"]))

    def f(name, step):
        return z[f"{name}__step{step}"]

    worst = {}
    for step in range(1, meta["nsteps"] + 1):
        prev = step - 1
        state = dict(
            T_mid=f("T_mid", prev), qv=f("qv", prev), qc=f("qc", prev),
            tke=f("tke", prev), u_wind=f("horiz_winds", prev)[:, 0, :],
            v_wind=f("horiz_winds", prev)[:, 1, :],
            cldfrac_liq=f("cldfrac_liq", prev),
            sgs_buoy_flux=f("sgs_buoy_flux", prev),
            eddy_diff_mom=f("eddy_diff_mom", prev))
        forcings = dict(
            p_mid=f("p_mid", prev), p_int=f("p_int", prev),
            pseudo_density=f("pseudo_density", prev), omega=f("omega", prev),
            phis=f("phis", prev), surf_sens_flux=f("surf_sens_flux", prev),
            surf_evap=f("surf_evap", prev),
            surf_mom_flux_x=f("surf_mom_flux", prev)[:, 0],
            surf_mom_flux_y=f("surf_mom_flux", prev)[:, 1])

        out = None
        for _ in range(N_SUBCYCLES):
            out = shoc_process_step(
                dt, npbl, cell_length,
                PARAMS["lambda_low"], PARAMS["lambda_high"],
                PARAMS["lambda_slope"], PARAMS["lambda_thresh"],
                PARAMS["thl2tune"], PARAMS["qw2tune"], PARAMS["qwthl2tune"],
                PARAMS["w2tune"], PARAMS["length_fac"],
                PARAMS["c_diag_3rd_mom"], PARAMS["ckh"], PARAMS["ckm"],
                False, False,
                state["T_mid"], forcings["p_mid"], forcings["p_int"],
                forcings["pseudo_density"], forcings["omega"],
                forcings["phis"], forcings["surf_sens_flux"],
                forcings["surf_evap"], forcings["surf_mom_flux_x"],
                forcings["surf_mom_flux_y"],
                state["qv"], state["qc"], state["tke"],
                state["u_wind"], state["v_wind"], state["cldfrac_liq"],
                state["sgs_buoy_flux"], state["eddy_diff_mom"])
            for k in PROGNOSTIC_MAP:
                state[k] = np.asarray(out[PROGNOSTIC_MAP[k]])
            state["u_wind"] = np.asarray(out["u_wind"])
            state["v_wind"] = np.asarray(out["v_wind"])

        # Compare against the captured post-step state
        def rel_err(mine, ref):
            scale = max(np.abs(ref).max(), 1e-30)
            return np.abs(mine - ref).max() / scale

        errs = {}
        for name in PROGNOSTIC_MAP:
            errs[name] = rel_err(state[name], f(name, step))
        uv_ref = f("horiz_winds", step)
        errs["u_wind"] = rel_err(state["u_wind"], uv_ref[:, 0, :])
        errs["v_wind"] = rel_err(state["v_wind"], uv_ref[:, 1, :])
        for name, key in DIAGNOSTIC_MAP.items():
            gname = f"{name}__step{step}"
            if gname not in z.files:
                continue
            if name == "inv_qc_relvar":
                # This field branches on shoc_ql2 != 0 exactly; points where
                # ql2 sits within ~1e-13 of zero flip the branch between the
                # C++ and the port. Require the flip fraction to be tiny and
                # everything else to match tightly.
                mine, ref = np.asarray(out[key]), z[gname]
                scale = max(np.abs(ref).max(), 1e-30)
                mismatch = np.abs(mine - ref) / scale > 1e-6
                assert mismatch.mean() < 1e-3, (
                    f"inv_qc_relvar branch-flip fraction too large: "
                    f"{mismatch.mean():.2e}")
                errs[name] = rel_err(mine[~mismatch], ref[~mismatch])
            else:
                errs[name] = rel_err(np.asarray(out[key]), z[gname])

        for k, v in errs.items():
            worst[k] = max(worst.get(k, 0.0), v)

    report = "\n".join(f"  {k:20s} max rel err = {v:.3e}"
                       for k, v in sorted(worst.items()))
    print(f"\nTier-1 SHOC golden replay ({meta['nsteps']} steps):\n{report}")

    bad = {k: v for k, v in worst.items() if v > 1e-6}
    assert not bad, f"fields exceeding 1e-6 relative error:\n{report}"
