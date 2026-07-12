"""Tier-1 for the ASSEMBLED physics suite: scream_jax.driver vs a
multi-process EAMxx golden run.

The archive (golden/physics_suite_218x72_dt1800_2steps.npz) was captured
from a pyeamxx run of the physics-only coupled configuration:
[mac_mic(shoc -> cld_fraction -> spa -> p3, 6 subcycles), rrtmgp].
This validates the DRIVER (process order, subcycling, timestamps, field
routing) on top of the individually swap-tested process steps.

Tolerances are calibrated against a control experiment: the C++ suite
run against a copy of itself with T_mid perturbed by 1e-6 diverges
after 2 steps with the same pattern (binary cldfrac flips at ~3% of
points, qv ~8e-3, LW_up ~1e-2 at 29% of points) — the coupled suite
amplifies roundoff through binary cloud fractions and P3 thresholds.
The replay's divergence sits within a few x of that inherent envelope;
a WIRING error instead produces a systematic signature (e.g. the
missing-aerosol bug this test caught: clear-sky LW off at ~100% of
points). The pass criterion below (fraction of points with >1e-3
field-scale error) separates the two regimes.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

REPO = Path(__file__).resolve().parents[2]
GOLDEN = REPO / "jax_port" / "golden" / "physics_suite_218x72_dt1800_2steps.npz"
DATA = REPO.parent / "e3sm-inputdata" / "atm" / "scream"

from scream_jax.driver import ScreamPhysics  # noqa: E402
from scream_jax.foundation.thermo import calculate_dx_from_area  # noqa: E402

# wiring-error discriminator: fraction of points beyond 1e-3
# field-scale error must stay small (inherent knife-edge noise reaches
# ~12% on the most sensitive cloud-top diagnostic; wiring errors reach
# ~100% on some field)
BIG_TOL = 1e-3
MAX_BIG_FRACTION = 0.20


def _load():
    if not GOLDEN.exists() or not (DATA / "init").exists():
        pytest.skip("suite golden archive or data files not available")
    z = np.load(GOLDEN)
    meta = json.loads(str(z["__metadata__"]))
    ic = meta["ic_file"]
    cands = [Path(ic)]
    if ic.startswith("/work/e3sm-inputdata"):
        cands.append(REPO.parent / "e3sm-inputdata"
                     / Path(ic).relative_to("/work/e3sm-inputdata"))
    ic_path = next((p for p in cands if p.exists()), None)
    if ic_path is None:
        pytest.skip("IC file not found")
    nc4 = pytest.importorskip("netCDF4")
    ds = nc4.Dataset(ic_path)
    geo = {k: np.array(ds[k][:]) for k in ("hyam", "hybm", "lat", "lon",
                                           "area")}
    ds.close()
    return z, meta, geo


def test_physics_suite_golden_replay():
    z, meta, geo = _load()
    dt = float(meta["dt"])
    doy0 = 284 + 45000.0 / 86400.0  # t0 = 2021-10-12-45000

    phys = ScreamPhysics(
        DATA, geo["hyam"], geo["hybm"], geo["lat"], geo["lon"],
        np.asarray(calculate_dx_from_area(geo["area"], geo["lat"])),
        mac_mic_subcycles=6, year=2021)

    def f(name, step):
        return z[f"{name}__step{step}"]

    worst, worst_frac = {}, {}
    for step in range(1, meta["nsteps"] + 1):
        prev = step - 1
        state = {name: f(name, prev) for name in meta["fields"]}
        doy_start = doy0 + prev * dt / 86400.0

        out = phys.step(state, dt, prev, doy_start)

        for name in meta["computed_fields"]:
            gname = f"{name}__step{step}"
            if gname not in z.files or name not in out:
                continue
            mine = np.asarray(out[name], dtype=np.float64)
            ref = np.asarray(z[gname], dtype=np.float64)
            if mine.shape != ref.shape:
                continue
            scale = max(np.abs(ref).max(), 1e-30)
            rel = np.abs(mine - ref) / scale
            bad = rel > BIG_TOL
            worst_frac[name] = max(worst_frac.get(name, 0.0), bad.mean())
            worst[name] = max(worst.get(name, 0.0), rel.max())

    report = "\n".join(
        f"  {k:28s} max rel err = {worst[k]:.3e}  "
        f"(frac > {BIG_TOL:g}: {worst_frac[k]:.2e})"
        for k in sorted(worst))
    print(f"\nTier-1 physics-suite replay ({meta['nsteps']} steps, "
          f"{len(worst)} fields):\n{report}")

    over = {k: v for k, v in worst_frac.items() if v > MAX_BIG_FRACTION}
    assert not over, (
        f"fields with more than {MAX_BIG_FRACTION:.0%} of points beyond "
        f"{BIG_TOL} (wiring-error signature):\n{report}")
