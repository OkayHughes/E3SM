"""Tier-1 golden replay + Tier-0 properties for the ZM convective
transport routines (eam_jax.zm_transport).

Golden archive: harness/gen_zm_transport_golden.py (run in the
scream-dev container). The mu/md/du/eu/ed/dp mass fluxes stored there
are physically consistent ZM plume structures: the Fortran
eam_zm_conv_f extension re-run on the zm_conv golden profiles (EAMv3
phys="default" params, config a; 18 of 42 columns triggered), plus
synthetic positive-definite tracers with contrasting vertical
structures, a 'dry'-type tracer, a negative-valued tracer, and
sheared/uniform momentum profiles. jt/mx/ideep are stored 1-based
(Fortran); the port is 0-based, hence the -1 shifts.

Tolerances (all measured):
- tracer dqdt: rtol 1e-12 (measured max 6.5e-13, 1-ulp cross-libm
  noise in the chat geometric log-average), atol 1e-24 (the
  |netflux| < max(flux)*1e-12 clip leaves denormal-scale entries;
  measured max abs diff 4.4e-25, physical tendencies are >= 1e-15).
- momentum icwu/icwd/pguall/pgdall: rtol 2e-12 (measured 2.4e-13; pg
  terms replay bitwise).
- momentum wind_tend/seten: rtol 1e-11, atol 1e-16. These are
  near-cancelling flux differences mu*(wiu-wi); measured rel 6.6e-13
  (wind_tend) / 1.9e-12 (seten) on the sheared case, and pure
  roundoff noise ~1e-18 on the uniform-wind case (covered by the
  atol, physical values are >= 1e-6).
"""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
GOLDEN = ROOT / "golden" / "zm_transport_golden.npz"

from eam_jax import zm_transport  # noqa: E402

TCFGS = ["t1", "t2", "t3"]
MCFGS = ["mom1", "mom2"]
MNAMES = ("wind_tend", "pguall", "pgdall", "icwu", "icwd", "seten")


@pytest.fixture(scope="module")
def gold():
    if not GOLDEN.exists():
        pytest.skip("golden archive missing (run harness in container)")
    return np.load(GOLDEN)


@pytest.fixture(scope="module")
def meta(gold):
    return json.loads(str(gold["__metadata__"]))


def _flux_args(gold):
    n = int(gold["lengath"])
    return dict(mu=gold["mu"], md=gold["md"], du=gold["du"],
                eu=gold["eu"], ed=gold["ed"], dp=gold["dp"],
                jt=gold["jt"] - 1, mx=gold["mx"] - 1,
                ideep=gold["ideep"][:n] - 1)


@pytest.fixture(scope="module")
def tresults(gold, meta):
    out = {}
    for cfg in TCFGS:
        c = meta["tracer_configs"][cfg]
        out[cfg] = np.asarray(zm_transport.zm_transport_tracer(
            gold["q"], gold["fracis"], np.array(c["doconvtran"]) != 0,
            gold["is_dry"] != 0, dpdry=gold["dpdry"], dt=c["dt"],
            zm_microp=bool(c["zm_microp"]), **_flux_args(gold)))
    return out


@pytest.fixture(scope="module")
def mresults(gold):
    out = {}
    for cfg in MCFGS:
        r = zm_transport.zm_transport_momentum(
            gold[f"wind_in_{cfg}"], dt=float(gold[f"dt_{cfg}"]),
            **_flux_args(gold))
        out[cfg] = {k: np.asarray(v) for k, v in r.items()}
    return out


# ---------------------------------------------------------------------------
# Tier-1 golden replay
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("cfg", TCFGS)
def test_replay_tracer(gold, tresults, cfg):
    np.testing.assert_allclose(tresults[cfg], gold[f"dqdt_{cfg}"],
                               rtol=1e-12, atol=1e-24)


@pytest.mark.parametrize("cfg", MCFGS)
@pytest.mark.parametrize("field", MNAMES)
def test_replay_momentum(gold, mresults, cfg, field):
    rtol, atol = (2e-12, 0.0) if field in ("icwu", "icwd", "pguall",
                                           "pgdall") \
        else (1e-11, 1e-16)  # near-cancelling flux differences
    np.testing.assert_allclose(mresults[cfg][field],
                               gold[f"{field}_{cfg}"],
                               rtol=rtol, atol=atol)


# ---------------------------------------------------------------------------
# Tier-0 properties
# ---------------------------------------------------------------------------
def test_untriggered_columns_exact_zero(gold, tresults, mresults):
    """Columns outside the gather set produce exactly zero tendencies;
    in-cloud winds equal the environment there."""
    n = int(gold["lengath"])
    ncol = gold["q"].shape[0]
    active = set((gold["ideep"][:n] - 1).tolist())
    idle = sorted(set(range(ncol)) - active)
    assert idle and active
    for cfg in TCFGS:
        assert np.all(tresults[cfg][idle] == 0.0)
    for cfg in MCFGS:
        r = mresults[cfg]
        for f in ("wind_tend", "pguall", "pgdall", "seten"):
            assert np.all(r[f][idle] == 0.0)
        np.testing.assert_array_equal(r["icwu"][idle],
                                      gold[f"wind_in_{cfg}"][idle])
        np.testing.assert_array_equal(r["icwd"][idle],
                                      gold[f"wind_in_{cfg}"][idle])


def test_tracer_column_mass_conservation(gold, meta, tresults):
    """sum_k dqdt*dptmp == 0 per gathered column to solver tolerance
    (the min(chat,const) limited fluxes telescope; the only
    non-conservative path is the zm_microp fixer, excluded here)."""
    n = int(gold["lengath"])
    gset = gold["ideep"][:n] - 1
    dp, dpdry = gold["dp"][:n], gold["dpdry"][:n]
    for cfg in ("t1", "t3"):
        doconv = meta["tracer_configs"][cfg]["doconvtran"]
        dqdt = tresults[cfg]
        for m in range(1, gold["q"].shape[2]):
            if not doconv[m]:
                continue
            dpm = dpdry if gold["is_dry"][m] else dp
            col = np.abs((dqdt[gset, :, m] * dpm).sum(axis=1))
            scale = max(np.abs(dqdt[gset, :, m] * dpm).max(), 1e-30)
            assert col.max() < 1e-10 * scale, (cfg, m)


def test_fixer_positivity(gold, meta, tresults):
    """With zm_microp the conservation fixer guarantees
    q + dt*dqdt >= 0 (to roundoff) for positive-definite tracers; the
    same tracer goes negative without it (the golden generator
    verified the fixer actually fires)."""
    dt = meta["tracer_configs"]["t2"]["dt"]
    q7 = gold["q"][:, :, 6]
    assert (q7 + dt * tresults["t1"][:, :, 6]).min() < 0.0
    # every positive-definite transported tracer stays nonnegative
    # after a t2 step (m=3 has negative INPUT values and is excluded)
    for m in (1, 2, 4, 5, 6):
        qn = gold["q"][:, :, m] + dt * tresults["t2"][:, :, m]
        assert qn.min() >= -1e-20, m
    # the uniform tracer never goes negative, so the fixer leaves its
    # (zero) tendency untouched
    np.testing.assert_array_equal(tresults["t1"][:, :, 4],
                                  tresults["t2"][:, :, 4])


def test_uniform_tracer_noop(gold, tresults):
    """A vertically uniform, fully insoluble tracer is transported
    without tendency (plume mass continuity): m=5 is constant 1e-6."""
    for cfg in TCFGS:
        assert np.abs(tresults[cfg][:, :, 4]).max() == 0.0


def test_water_vapor_and_disabled_slices_zero(gold, meta, tresults):
    """m=0 (water vapor) is never transported here; doconvtran=False
    slices stay exactly zero."""
    for cfg in TCFGS:
        doconv = meta["tracer_configs"][cfg]["doconvtran"]
        assert np.all(tresults[cfg][:, :, 0] == 0.0)
        for m, on in enumerate(doconv):
            if not on:
                assert np.all(tresults[cfg][:, :, m] == 0.0), (cfg, m)


def test_momentum_column_conservation(gold, mresults):
    """sum_k wind_tend*dp == 0 per gathered column (the momentum flux
    divergence telescopes to the zero cloud-top flux)."""
    n = int(gold["lengath"])
    gset = gold["ideep"][:n] - 1
    dp = gold["dp"][:n]
    for cfg in MCFGS:
        wt = mresults[cfg]["wind_tend"]
        for m in (0, 1):
            col = np.abs((wt[gset, :, m] * dp).sum(axis=1))
            assert col.max() < 1e-13, (cfg, m)


def test_uniform_wind_noop(gold, mresults):
    """Uniform winds (mom2): no shear -> no pressure-gradient force,
    in-cloud winds equal the environment, zero tendency and zero KE
    dissipation to roundoff."""
    r = mresults["mom2"]
    assert np.abs(r["wind_tend"]).max() < 1e-16
    assert np.abs(r["seten"]).max() < 1e-15
    assert np.abs(r["pguall"]).max() == 0.0
    np.testing.assert_allclose(r["icwu"], gold["wind_in_mom2"],
                               rtol=1e-13)


def test_ke_dissipation_energy_closure(gold, mresults):
    """B&B03 energy closure: seten = ketend_cons - ketend, where the
    conservative part telescopes to the (zero) cloud-top flux, so the
    column integral of seten must equal minus the column-integrated
    KE change of the updated winds: sum(seten*dp) == -sum(0.5*
    (wf^2 - w0^2)/dt * dp) with wf = w0 + dt*wind_tend (identical to
    the scheme's windf up to roundoff). Measured residual 8.9e-16 on
    column integrals of order 2; asserted at atol 1e-12."""
    n = int(gold["lengath"])
    gset = gold["ideep"][:n] - 1
    dp = gold["dp"][:n]
    for cfg in MCFGS:
        r = mresults[cfg]
        w0 = gold[f"wind_in_{cfg}"][gset]
        dt = float(gold[f"dt_{cfg}"])
        col = (r["seten"][gset] * dp).sum(axis=1)
        wf = w0 + dt * r["wind_tend"][gset]
        ke_change = (((wf ** 2).sum(axis=2) - (w0 ** 2).sum(axis=2))
                     * 0.5 / dt * dp).sum(axis=1)
        np.testing.assert_allclose(col, -ke_change,
                                   rtol=1e-12, atol=1e-12)
