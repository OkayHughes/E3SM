"""P3 constants, verbatim from micro_p3_utils.F90 (eam variant) and the
physconst values fed to micro_p3_utils_init by micro_p3_interface.F90
(physconst derives them from share/util/shr_const_mod.F90 — see
eam_jax/constants.py, the shared source of truth for these).

Single-precision literal promotions in the Fortran are reproduced with
np.float32 round-trips and commented at each site.
"""

import numpy as np

from ..constants import (  # physconst-derived, verbatim shr_const
    SHR_CONST_PI,
    SHR_CONST_RGAS,
)

# micro_p3_utils_init arguments (physconst)
pi = SHR_CONST_PI
cp = 1.00464e3                     # cpair
inv_cp = 1.0 / cp
g = 9.80616                        # gravit
rd = SHR_CONST_RGAS / 28.966       # rair
rv = SHR_CONST_RGAS / 18.016       # rh2o
ep_2 = 18.016 / 28.966             # mwh2o/mwdry
rho_h2o = 1.000e3                  # rhoh2o
cpw = 4.188e3                      # cpliq
T_zerodegc = 273.15                # tmelt
T_homogfrz = T_zerodegc - 40.0
T_icenuc = T_zerodegc - 15.0
T_rainfrz = T_zerodegc - 4.0
latvap = 2.501e6
latice = 3.337e5
latsub = latvap + latice

rho_1000mb = 100000.0 / (rd * T_zerodegc)
rho_600mb = 60000.0 / (rd * 253.15)
inv_rho_h2o = 1.0 / rho_h2o
dropmass = 5.2e-7
inv_dropmass = 1.0 / dropmass

# module parameters (micro_p3_utils.F90)
qsmall = 1.0e-14
nsmall = 1.0e-16
thrd = 1.0 / 3.0
sxth = 1.0 / 6.0
piov3 = pi * thrd
piov6 = pi * sxth
max_total_ni = 500.0e3             # EAM value (SCREAM uses 740e3)
kc = 9.44e9
kr = 5.78e3
ar = 841.99667
br = 0.8
f1r = 0.78
f2r = 0.32
ecr = 1.0
rho_rimeMin = 50.0
rho_rimeMax = 900.0
inv_rho_rimeMax = 1.0 / rho_rimeMax
bimm = 2.0
aimm = 0.65
rin = 0.1e-6
mi0 = 4.0 * piov3 * 900.0 * 1.0e-18
eci = 0.5
eri = 1.0
bcn = 2.0
dbrk = 600.0e-6
nmltratio = 1.0
clbfact_dep = 1.0
clbfact_sub = 1.0
mu_r_constant = 0.0                # EAM value (SCREAM uses 1.0)
lookup_table_1a_dum1_c = 4.135985029041767  # 1/(0.1*log10(261.7))
mincld = 0.0001
rho_h2os = 917.0
iparam = 3                         # Khairoutdinov & Kogan 2000

min_mean_mass_liq = 1.0e-20
min_mean_mass_ice = 1.0e-20
min_cld_frac = 1.0e-20
# single-precision literals in micro_p3_utils.F90 (5.1E-3, 1.0E-2)
incloud_limit = float(np.float64(np.float32(5.1e-3)))
precip_limit = float(np.float64(np.float32(1.0e-2)))

cons1 = piov6 * rho_h2o
cons2 = 4.0 * piov3 * rho_h2o
cons3 = 1.0 / cons2                # embryonic size applied at call site
cons4 = 1.0 / (dbrk ** 3 * pi * rho_h2o)
cons5 = piov6 * bimm
cons6 = piov6 ** 2 * rho_h2o * bimm
cons7 = 4.0 * piov3 * rho_h2o * 1.0e-18

# droplet spectral shape parameter table (iparam=1 only; kept verbatim)
dnu = np.array([0.000, -0.557, -0.430, -0.307, -0.186, -0.067,
                -0.050, -0.167, -0.282, -0.397, -0.512, -0.626,
                -0.739, -0.853, -0.966, -0.966])

# ice lookup table dimensions
isize, densize, rimsize, rcollsize = 50, 4, 5, 30   # noqa: E501 (isize, rimsize order below)
isize = 50
densize = 5
rimsize = 4
rcollsize = 30
ice_table_size = 12
collect_table_size = 2

# ice_complete_melting threshold: 273.15 is a DEFAULT-REAL (single)
# literal in the Fortran (t_snow_melt = 273.15 + 2.0_rtype)
t_snow_melt = float(np.float64(np.float32(273.15)) + 2.0)

# sedimentation substep tolerance (do while (dt_left .gt. 1.e-4))
dt_left_tol = 1.0e-4
