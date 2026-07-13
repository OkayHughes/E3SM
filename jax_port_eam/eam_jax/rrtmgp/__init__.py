"""JAX port of the EAM (E3SMv3) RRTMGP radiation layer.

PORT_NOTES — module provenance (see each module's docstring for the
line-level notes):

  coefficients.py  COPIED from ../jax_port/scream_jax/rrtmgp (identical
                   algorithm; the C++ EAMxx loader is a direct port of
                   the Fortran mo_load_coefficients + gas-optics load).
  gas_optics.py    COPIED (Fortran kernels algorithmically identical to
                   the EAMxx C++; constants match mo_rrtmgp_constants
                   defaults, which EAM never overrides).
  rte.py           COPIED + net-flux outputs added (EAM's fluxes_t
                   carries flux_net = dn - up).
  optical_props.py COPIED + delta_scale adapted to the FORTRAN kernel
                   (no tau gate, max(eps, .) guarded divisions).
  mcica.py         FRESH (EAM uses the share/RandNum KISS generator
                   seeded from bottom-layer pressures; EAMxx uses JSF64
                   seeded per cell — completely different streams).
  cloud_optics.py  FRESH (EAM default gammadist liquid + mitchell ice
                   tables with the CAM lininterp; EAMxx uses the RRTMGP
                   cloud-optics LUT).
  driver.py        FRESH (radiation.F90 radiation_tend sequencing +
                   the f90 rrtmgp_interface run wrappers).

Level ordering: level 0 = model top everywhere (EAM and this port);
the radiation grid adds one level above the model top.
"""
