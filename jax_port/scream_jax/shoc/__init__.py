"""SHOC turbulence/macrophysics (shoc process).

Transcribed kernel-by-kernel from
components/eamxx/src/physics/shoc/impl/*_impl.hpp (master header
shoc_functions.hpp). Modules are grouped by subsystem following the
call-tree in jax_port/PORTING_PLAN.md section 4.3; each function maps 1:1
to a C++ kernel and names match the C++ (minus the shoc_ prefix where the
module name already carries it).

Conventions (identical to EAMxx/C++):
- k=0 is the model top. Midpoint ("zt") arrays have nlev entries,
  interface ("zi") arrays nlevi = nlev+1. The level axis is the LAST axis;
  leading axes are batch (columns).
- Kernels that update only a band of levels (e.g. interior interfaces)
  take the current array and return a functionally-updated copy.
"""

from . import constants
from .assumed_pdf import shoc_assumed_pdf
from .main import shoc_init, shoc_main
from .energy import (
    shoc_energy_fixer,
    shoc_energy_integrals,
    update_host_dse,
)
from .grid import shoc_grid, dp_inverse, compute_tmpi
from .interp import linear_interp
from .length import (
    check_length_scale_shoc_length,
    compute_brunt_shoc_length,
    compute_l_inf_shoc_length,
    compute_shoc_mix_shoc_length,
    shoc_length,
)
from .pblintd import (
    pblintd,
    pblintd_check_pblh,
    pblintd_height,
    pblintd_surf_temp,
    shoc_pblintd_cldcheck,
    shoc_pblintd_init_pot,
)
from .second_moments import (
    calc_shoc_vertflux,
    calc_shoc_varorcovar,
    diag_second_moments,
    diag_second_shoc_moments,
    shoc_diag_second_moments_lbycond,
    shoc_diag_second_moments_srf,
)
from .solver import (
    update_prognostics_implicit,
    vd_shoc_decomp,
    vd_shoc_solve,
)
from .surface import shoc_diag_obklen
from .thermo import compute_shoc_vapor, compute_shoc_temperature
from .third_moments import (
    clipping_diag_third_shoc_moments,
    compute_diag_third_shoc_moment,
    diag_third_shoc_moments,
)
from .tke import (
    adv_sgs_tke,
    check_tke,
    compute_shr_prod,
    eddy_diffusivities,
    integ_column_stability,
    isotropic_ts,
    shoc_tke,
)
