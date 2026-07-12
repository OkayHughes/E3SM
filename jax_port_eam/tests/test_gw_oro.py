"""Tier-1 golden replay + Tier-0 properties for the orographic
gravity-wave spine (eam_jax.gw)."""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
GOLDEN = ROOT / "golden" / "gw_oro_golden.npz"

from eam_jax import gw  # noqa: E402


@pytest.fixture(scope="module")
def gold():
    if not GOLDEN.exists():
        pytest.skip("golden archive missing (run harness in container)")
    return np.load(GOLDEN)


def params_from(gold, oro_only):
    meta = json.loads(str(gold["__metadata__"]))
    p = meta["params"]
    return gw.make_gw_params(gold["alpha"], p["kbotbg"], meta["nlev"],
                             fcrit2=p["fcrit2"], kwv=p["kwv"],
                             gravit=p["gravit"], rair=p["rair"],
                             ktop=p["ktop"], orographic_only=oro_only)


@pytest.mark.parametrize("oro", [0, 1])
def test_gw_prof(gold, oro):
    params = params_from(gold, bool(oro))
    rhoi, ti, nm, ni = gw.gw_prof(1004.64, gold["t"], gold["pmid"],
                                  gold["pint"], params)
    s = f"_oro{oro}"
    np.testing.assert_allclose(np.asarray(rhoi), gold["rhoi" + s],
                               rtol=1e-13)
    np.testing.assert_allclose(np.asarray(ti), gold["ti" + s], rtol=1e-13)
    np.testing.assert_allclose(np.asarray(nm), gold["nm" + s], rtol=1e-13)
    np.testing.assert_allclose(np.asarray(ni), gold["ni" + s], rtol=1e-13)


@pytest.mark.parametrize("oro", [0, 1])
def test_gw_oro_src(gold, oro):
    params = params_from(gold, bool(oro))
    s = f"_oro{oro}"
    src, tend, tau, ubm, ubi, xv, yv = gw.gw_oro_src(
        gold["u"], gold["v"], gold["t"], gold["sgh"], gold["pmid"],
        gold["pint"], gold["dpm"], gold["zm"], gold["nm" + s], params)
    np.testing.assert_array_equal(np.asarray(src), gold["src_level" + s])
    np.testing.assert_array_equal(np.asarray(tend),
                                  gold["tend_level" + s])
    np.testing.assert_allclose(np.asarray(tau), gold["tau0" + s],
                               rtol=1e-12, atol=1e-300)
    np.testing.assert_allclose(np.asarray(ubm), gold["ubm" + s],
                               rtol=1e-12, atol=1e-14)
    np.testing.assert_allclose(np.asarray(ubi), gold["ubi" + s],
                               rtol=1e-12, atol=1e-14)
    np.testing.assert_allclose(np.asarray(xv), gold["xv" + s], rtol=1e-13)
    np.testing.assert_allclose(np.asarray(yv), gold["yv" + s], rtol=1e-13)


@pytest.mark.parametrize("oro", [0, 1])
def test_gw_drag_prof(gold, oro):
    meta = json.loads(str(gold["__metadata__"]))
    params = params_from(gold, bool(oro))
    s = f"_oro{oro}"
    tau, utgw, vtgw = gw.gw_drag_prof_oro(
        gold["src_level" + s], gold["tend_level" + s], meta["dt"],
        gold["lat"], gold["t"], gold["ti" + s], gold["pmid"],
        gold["pint"], gold["dpm"], gold["rdpm"], gold["piln"],
        gold["rhoi" + s], gold["nm" + s], gold["ni" + s],
        gold["ubm" + s], gold["ubi" + s], gold["xv" + s],
        gold["yv" + s], meta["params"]["effgw_oro"], gold["tau0" + s],
        params)
    np.testing.assert_allclose(np.asarray(tau), gold["tau" + s],
                               rtol=1e-12, atol=1e-300)
    np.testing.assert_allclose(np.asarray(utgw), gold["utgw" + s],
                               rtol=1e-12, atol=1e-20)
    np.testing.assert_allclose(np.asarray(vtgw), gold["vtgw" + s],
                               rtol=1e-12, atol=1e-20)


# ---------------------------------------------------------------------------
# Tier-0 properties
# ---------------------------------------------------------------------------
def test_ocean_columns_zero_stress(gold):
    # sgh = 0 columns (0..3) must produce zero stress and tendency
    params = params_from(gold, True)
    s = "_oro1"
    assert gold["tau0" + s][:4].max() == 0.0
    np.testing.assert_array_equal(gold["utgw" + s][:4], 0.0)


def test_tendency_decelerates_source_wind(gold):
    # oro c=0 drag opposes the projected source wind: gwut has sign of
    # (0 - ubm); check the projected tendency is anti-correlated with ubm
    s = "_oro1"
    ub_tend = (gold["utgw" + s] * gold["xv" + s][:, None]
               + gold["vtgw" + s] * gold["yv" + s][:, None])
    mask = np.abs(ub_tend) > 0
    assert (ub_tend[mask] * gold["ubm" + s][mask] <= 0).all()


def test_tendency_limit(gold):
    # |projected tendency| <= tndmax * effgw (oro path)
    s = "_oro1"
    ub_tend = np.hypot(gold["utgw" + s], gold["vtgw" + s])
    assert ub_tend.max() <= 500.0 / 86400.0 * 0.375 * (1 + 1e-12)
