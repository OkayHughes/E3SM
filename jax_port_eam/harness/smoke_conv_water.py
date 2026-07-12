#!/usr/bin/env python3
"""Container smoke test for the conv_water f2py extension: float64
wrapper check plus a hand-computed mode-1 case and branch sanity."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "fmod"))
import eam_conv_water_f  # noqa: E402

dcw = eam_conv_water_f.conv_water_driver.drv_conv_water
assert "array('d')" in dcw.__doc__, "wrapper is not float64!"
print("conv_water wrapper float64 OK")

ncol, nlev = 4, 72
F = lambda a: np.asarray(a, dtype=np.float64, order="F")  # noqa: E731
t = F(np.full((ncol, nlev), 280.0))
pdel = F(np.full((ncol, nlev), 1000.0))
ql = F(np.full((ncol, nlev), 2e-5))
qi = F(np.full((ncol, nlev), 1e-5))
shw = F(np.full((ncol, nlev), 1e-3))
dpw = F(np.full((ncol, nlev), 2e-3))
dpi = F(np.full((ncol, nlev), 5e-4))
fice = F(np.full((ncol, nlev), 0.3))
shf = F(np.full((ncol, nlev), 0.05))
dpf = F(np.full((ncol, nlev), 0.10))
ast = F(np.full((ncol, nlev), 0.4))
rei = F(np.full((ncol, nlev), 25.0))

tl, ti, scl, sci = dcw(1, 0, 0, 0, t, pdel, ql, qi, shw, dpw, dpi,
                       fice, shf, dpf, ast, rei)
# hand check: mode 1, zm_microp=F
wrk1 = 1e-5 / (3e-5 + 1e-36)
cu0 = 0.15
cu_ic = (0.05 * 1e-3 + 0.10 * 2e-3) / cu0
ls_ic = 3e-5 / 0.4
tot0 = 0.55
tot_ic = (0.4 * ls_ic + cu0 * cu_ic) / tot0
exp_ti = tot0 * tot_ic * wrk1
exp_tl = tot0 * tot_ic * (1 - wrk1)
print("mode1 rel err liq/ice:", abs(tl[0, 0] - exp_tl) / exp_tl,
      abs(ti[0, 0] - exp_ti) / exp_ti)
assert abs(tl[0, 0] - exp_tl) / exp_tl < 1e-14
assert abs(ti[0, 0] - exp_ti) / exp_ti < 1e-14
assert abs(scl[0, 0] + sci[0, 0] - 1e-3 * 0.05) < 1e-18

# NaN fice guards the COSP outputs
fice2 = fice.copy(order="F")
fice2[1, :] = np.nan
_, _, scl, sci = dcw(1, 0, 0, 0, t, pdel, ql, qi, shw, dpw, dpi,
                     fice2, shf, dpf, ast, rei)
assert np.all(scl[1, :] == 0) and np.all(sci[1, :] == 0)

# zm_microp branch ignores conv_water_mode
o1 = dcw(1, 1, 0, 0, t, pdel, ql, qi, shw, dpw, dpi, fice, shf, dpf,
         ast, rei)
o2 = dcw(2, 1, 0, 0, t, pdel, ql, qi, shw, dpw, dpi, fice, shf, dpf,
         ast, rei)
for a, b in zip(o1, o2):
    assert np.array_equal(a, b)

# no clouds anywhere -> all zeros
z = F(np.zeros((ncol, nlev)))
o = dcw(1, 0, 0, 0, t, pdel, z, z, z, z, z, z, z, z, z, rei)
assert all(np.all(x == 0.0) for x in o)
print("conv_water smoke OK; microp totg_liq[0,0] =", o1[0][0, 0])
