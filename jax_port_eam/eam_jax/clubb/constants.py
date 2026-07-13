"""CLUBB model constants, verbatim from
components/eam/src/physics/clubb/constants_clubb.F90 under the EAM
build defines (-DCLUBB_CAM: physical constants come from the real
shr_const_mod, mirrored in eam_jax/constants.py)."""

from .. import constants as shr

# CLUBB_CAM physical constants (constants_clubb.F90 lines 183-217)
CP = shr.SHR_CONST_CPDAIR
LV = shr.SHR_CONST_LATVAP
LF = shr.SHR_CONST_LATICE
LS = shr.SHR_CONST_LATSUB
RD = shr.SHR_CONST_RDAIR
RV = shr.SHR_CONST_RGAS / shr.SHR_CONST_MWWV
T_FREEZE_K = shr.SHR_CONST_TKFRZ
EP = shr.SHR_CONST_MWWV / shr.SHR_CONST_MWDAIR   # 0.622
EP1 = (1.0 - EP) / EP                            # 0.61
EP2 = 1.0 / EP                                   # 1.61
KAPPA = (shr.SHR_CONST_RGAS / shr.SHR_CONST_MWDAIR) / shr.SHR_CONST_CPDAIR
GRAV = shr.SHR_CONST_G
P0 = 1.0e5
VONK = shr.SHR_CONST_KARMAN
RHO_LW = shr.SHR_CONST_RHOFW

# Numerical constants (verbatim literals)
SQRT_2PI = 2.5066282746310005024
SQRT_2 = 1.4142135623730950488

# Tolerances (constants_clubb.F90 lines 302-330)
W_TOL = 2.0e-2          # [m/s]
THL_TOL = 1.0e-2        # [K]
RT_TOL = 1.0e-8         # [kg/kg]
CHI_TOL = 1.0e-8        # [kg/kg]
ETA_TOL = CHI_TOL
W_TOL_SQD = W_TOL ** 2
RC_TOL = 1.0e-6         # [kg/kg]

EPS = 1.0e-10           # smallest divide-guard
MAX_NUM_STDEVS = 5.0    # range of stdevs for statistical significance
ZERO_THRESHOLD = 0.0
MAX_MAG_CORRELATION = 0.99

# Fortran epsilon(1.0_core_rknd) for double precision (clip_rcm uses
# rtm - epsilon(rtm)); 2**-52.
EPSILON_R8 = 2.0 ** -52
