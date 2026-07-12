! Stub physconst exposing the subset of constants the ported schemes
! use, derived VERBATIM from share/util/shr_const_mod.F90 (which is
! compiled from the real source) exactly as the real
! eam/src/utils/physconst.F90 derives them.
module physconst
  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_const_mod, only: &
      shr_const_mwdair, shr_const_mwwv, shr_const_rgas, &
      shr_const_rdair, shr_const_rwv, shr_const_zvir, &
      shr_const_cpdair, shr_const_cpwv, shr_const_latvap, &
      shr_const_latice, shr_const_latsub, shr_const_tkfrz, &
      shr_const_tktrip, shr_const_g, shr_const_stebol, &
      shr_const_karman, shr_const_avogad, shr_const_boltz, &
      shr_const_rhofw, shr_const_pstd, shr_const_pi, &
      shr_const_cpfw
  implicit none

  real(r8), public, parameter :: avogad = shr_const_avogad
  real(r8), public, parameter :: boltz = shr_const_boltz
  real(r8), public, parameter :: cpair = shr_const_cpdair
  real(r8), public, parameter :: cpliq = shr_const_cpfw
  real(r8), public, parameter :: cpwv = shr_const_cpwv
  real(r8), public, parameter :: epsilo = shr_const_mwwv/shr_const_mwdair
  real(r8), public, parameter :: gravit = shr_const_g
  real(r8), public, parameter :: h2otrip = shr_const_tktrip
  real(r8), public, parameter :: karman = shr_const_karman
  real(r8), public, parameter :: latice = shr_const_latice
  real(r8), public, parameter :: latvap = shr_const_latvap
  real(r8), public, parameter :: latsub = shr_const_latsub
  real(r8), public, parameter :: mwdry = shr_const_mwdair
  real(r8), public, parameter :: mwh2o = shr_const_mwwv
  real(r8), public, parameter :: pi = shr_const_pi
  real(r8), public, parameter :: rair = shr_const_rdair
  real(r8), public, parameter :: rga = 1._r8/shr_const_g
  real(r8), public, parameter :: rgas = shr_const_rgas
  real(r8), public, parameter :: rh2o = shr_const_rwv
  real(r8), public, parameter :: rhoh2o = shr_const_rhofw
  real(r8), public, parameter :: stebol = shr_const_stebol
  real(r8), public, parameter :: tmelt = shr_const_tkfrz
  real(r8), public, parameter :: zvir = shr_const_zvir
  real(r8), public, parameter :: pstd = shr_const_pstd
  real(r8), public, parameter :: cappa = shr_const_rdair/shr_const_cpdair
  real(r8), public, parameter :: cpvir = shr_const_cpwv/shr_const_cpdair - 1._r8
end module physconst
