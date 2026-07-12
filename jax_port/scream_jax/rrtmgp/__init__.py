"""RRTMGP radiation (rrtmgp process).

Transcribed from the Kokkos-templated RTE+RRTMGP implementation used by
EAMxx (components/eam/src/physics/rrtmgp/external/cpp, classes with a K
suffix) and the EAMxx interface layer
(components/eamxx/src/physics/rrtmgp/eamxx_rrtmgp_interface.hpp).

Conventions:
  - arrays are (col, lay[, gpt/bnd]) with k=0 at the model top unless
    noted; tables keep the C++ index order (e.g. kmajor(gpt, eta,
    press, temp)), obtained from the netCDF files by transposing the
    C-order data (the C++ conv::SimpleNetCDF reader reverses dims);
  - ALL integer index data read from files is shifted to 0-based, as in
    conv::SimpleNetCDF::read (it subtracts 1 from every int array);
  - band2gpt / minor_limits_gpt hold 0-based INCLUSIVE ranges.
"""
