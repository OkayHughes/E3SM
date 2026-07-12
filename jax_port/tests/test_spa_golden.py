"""Tier-1: scream_jax SPA vs the EAMxx golden archive.

Replays golden/spa_218x72_dt1800_3steps.npz: for each step, computes the
prescribed-aerosol fields at the end-of-step timestamp and compares
against the captured state.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

REPO = Path(__file__).resolve().parents[2]
GOLDEN = REPO / "jax_port" / "golden" / "spa_218x72_dt1800_3steps.npz"
SPA_FILE = REPO.parent / "e3sm-inputdata" / "atm" / "scream" / "init" / \
    "spa_file_unified_and_complete_ne2np4L72_20231222.nc"

from scream_jax.spa.process import (FIELD_MAP, load_spa_data,  # noqa: E402
                                    spa_process_step)


def test_spa_golden_replay():
    if not GOLDEN.exists() or not SPA_FILE.exists():
        pytest.skip("spa golden archive or data file not available")
    z = np.load(GOLDEN)
    meta = json.loads(str(z["__metadata__"]))
    dt = float(meta["dt"])
    # t0 = 2021-10-12-45000 -> 0-based fractional day-of-year
    doy0 = 284 + 45000.0 / 86400.0

    data = load_spa_data(str(SPA_FILE))

    worst = {}
    for step in range(1, meta["nsteps"] + 1):
        doy_end = doy0 + step * dt / 86400.0
        out = spa_process_step(data, doy_end, z[f"p_mid__step{step - 1}"])
        for name in FIELD_MAP:
            ref = z[f"{name}__step{step}"]
            mine = np.asarray(out[name])
            scale = max(np.abs(ref).max(), 1e-30)
            err = np.abs(mine - ref).max() / scale
            worst[name] = max(worst.get(name, 0.0), err)

    report = "\n".join(f"  {k:14s} max rel err = {v:.3e}"
                       for k, v in sorted(worst.items()))
    print(f"\nTier-1 SPA golden replay ({meta['nsteps']} steps):\n{report}")
    bad = {k: v for k, v in worst.items() if v > 1e-10}
    assert not bad, f"fields exceeding 1e-10 relative error:\n{report}"
