"""Physical and numerical constants used across the SCREAM physics suite.

Transcribed verbatim from
components/eamxx/src/share/physics/physics_constants.hpp
(struct scream::physics::Constants<Scalar>, double-precision instantiation).

Names and grouping follow the C++ source so a reader can diff this file
against it line by line. EKAT's compile-time unit annotations (quantity_t)
carry no runtime value and are dropped; only the SI magnitudes remain.
Values must NOT be "improved" (e.g. more digits of a physical constant):
agreement with EAMxx, not with CODATA, is the correctness criterion.
"""

import numpy as np

# --- Commonly used numeric constants -----------------------------------------
ZERO = 0.0
ONE = 1.0
THIRD = 1.0 / 3.0
SXTH = 1.0 / 6.0
Pi = 3.14159265358979323
PIOV3 = Pi * THIRD
PIOV6 = Pi * SXTH
macheps = float(np.finfo(np.float64).eps)

# --- Physical constants -------------------------------------------------------
Cpair = 1004.64          # heat capacity of dry air at constant pressure [J/kg/K]
CP = Cpair
INV_CP = ONE / CP
Rair = 287.042           # gas constant for dry air [J/kg/K]
RD = Rair
RH2O = 461.505           # water vapor gas constant [J/kg/K]
RV = RH2O

RHO_H2O = 1000.0         # liquid water density [kg/m3]
INV_RHO_H2O = 1.0 / RHO_H2O
RHOW = RHO_H2O
INV_RHOW = 1.0 / RHOW
RhoIce = 917.0           # ice density at 0 C, Wallace & Hobbs 1977 [kg/m3]

MWH2O = 18.016           # water molar mass [g/mol]
MWWV = MWH2O
MWdry = 28.966           # dry air molar mass [g/mol]
ep_2 = MWH2O / MWdry     # dimensionless
o2mmr = 0.23143          # O2 mass mixing ratio [-]

gravit = 9.80616         # gravitational acceleration [m/s2]

LatVap = 2501000.0       # latent heat of vaporization [J/kg]
LatIce = 333700.0        # latent heat of fusion [J/kg]
CpLiq = 4188.0           # heat capacity of liquid water [J/kg/K]

Tmelt = 273.15           # melting point [K]
T_zerodegc = Tmelt
T_homogfrz = Tmelt - 40.0
T_rainfrz = Tmelt - 4.0

RHO_RIMEMIN = 50.0       # min limit for rime density [kg/m3]
RHO_RIMEMAX = 900.0      # max limit for rime density [kg/m3]
INV_RHO_RIMEMAX = 1.0 / RHO_RIMEMAX
BIMM = 2.0
CONS1 = PIOV6 * RHOW
CONS2 = 4.0 * PIOV3 * RHOW
CONS3 = 1.0 / (CONS2 * 1.562500000000000e-14)  # 1/(CONS2*pow(25.e-6,3))
CONS5 = PIOV6 * BIMM
CONS6 = PIOV6 * PIOV6 * RHOW * BIMM
CONS7 = 4.0 * PIOV3 * RHOW * 1.0e-18
QSMALL = 1.0e-14
QTENDSMALL = 1.0e-20
BSMALL = 1.0e-15
NSMALL = 1.0e-16
P0 = 100000.0            # reference pressure [Pa]
RHOSUR = P0 / (RD * Tmelt)
rhosui = 60000.0 / (RD * 253.15)
RHO_1000MB = P0 / (RD * Tmelt)
RHO_600MB = 60000.0 / (RD * 253.15)
dt_left_tol = 1.0e-4
bcn = 2.0
dropmass = 5.2e-7
NCCNST = 200.0e6
incloud_limit = 5.1e-3
precip_limit = 1.0e-2
Karman = 0.4
Avogad = 6.02214e26      # [1/mol] (note: per kmol-style convention, as in EAMxx)
Boltz = 1.38065e-23      # [J/K]
Rgas = Avogad * Boltz
RWV = Rgas / MWWV
ZVIR = (RWV / Rair) - 1.0
f1r = 0.78
f2r = 0.32
nmltratio = 1.0          # ratio of rain number produced to ice number lost in melting
basetemp = 300.0
r_earth = 6.376e6        # Earth radius [m]
stebol = 5.670374419e-8  # Stefan-Boltzmann constant [W/m2/K4]
omega = 7.292e-5         # Earth's rotation rate [rad/s]

# --- Table dimension constants (P3 rain tables) --------------------------------
VTABLE_DIM0 = 300
VTABLE_DIM1 = 10
MU_R_TABLE_DIM = 150

# --- Turbulent Mountain Stress constants ---------------------------------------
orocnst = 1.0            # std deviation -> height conversion [-]
z0fac = 0.075            # z_0 from orographic standard deviation [-]

# Warm-rain parameterization switch: 1 = Seifert-Beheng 2001, 2 = Beheng 1994,
# 3 = Khairoutdinov-Kogan 2000
IPARAM = 3

# --- WGS84 ellipsoid coefficients (area -> length conversion) ------------------
earth_ellipsoid1 = 111132.92  # meters per degree longitude at equator
earth_ellipsoid2 = 559.82
earth_ellipsoid3 = 1.175

# --- Gas molecular weights (Constants::get_gas_mol_weight) ---------------------
_GAS_MOL_WEIGHTS = {
    "h2o": MWH2O,
    "co2": 44.0095,
    "o3": 47.9982,
    "n2o": 44.0128,
    "co": 28.0101,
    "ch4": 16.04246,
    "o2": 31.998,
    "n2": 28.0134,
    "cfc11": 136.0,
    "cfc12": 120.0,
}


def get_gas_mol_weight(gas_name: str) -> float:
    """Molecular weight [g/mol] of a radiatively active gas (case-insensitive)."""
    try:
        return _GAS_MOL_WEIGHTS[gas_name.lower()]
    except KeyError:
        raise ValueError(f"Unknown gas name: {gas_name!r}") from None
