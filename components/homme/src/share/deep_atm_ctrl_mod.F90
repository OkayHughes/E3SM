module deep_atm_mod
use physical_constants, only : rearth, gravit=>g, omega
use kinds, only: real_kind
use dimensions_mod, only: np
implicit none
private
save

public :: r_hat_from_phi, &
          z_from_phi, &
          g_from_phi, &
          g_from_z, &
          quasi_hydrostatic_terms
public :: atm_is_deep, gravity_is_const
#ifdef HOMMEDA
logical, parameter :: atm_is_deep=.true.
logical, parameter :: gravity_is_const=.false.
#else
logical, parameter :: atm_is_deep=.false.
logical, parameter :: gravity_is_const=.true.
#endif
interface r_hat_from_phi
         module procedure r_hat_from_phi_rank_ijk
         module procedure r_hat_from_phi_rank_ij 
        module procedure r_hat_from_phi_rank_k
         module procedure r_hat_from_phi_no_rank
end interface
interface z_from_phi
         module procedure z_from_phi_rank_ijk
         module procedure z_from_phi_rank_ij
         module procedure z_from_phi_rank_k
         module procedure z_from_phi_no_rank
end interface
interface g_from_phi
         module procedure g_from_phi_rank_ijk
         module procedure g_from_phi_rank_ij
         module procedure g_from_phi_rank_k
         module procedure g_from_phi_no_rank
end interface
interface g_from_z
         module procedure g_from_z_rank_ijk
         module procedure g_from_z_rank_ij
         module procedure g_from_z_rank_k
         module procedure g_from_z_no_rank
end interface




contains
function r_hat_from_phi_rank_ijk(phi, k) result(r_hat)
  integer,              intent(in) :: k
  real(kind=real_kind), intent(in) :: phi(np,np,k)
  real(kind=real_kind)           :: r_hat(np,np,k)
  r_hat = elementwise_ijk(phi, k, "r_hat_from_phi_scalar")
end function r_hat_from_phi_rank_ijk
function r_hat_from_phi_rank_ij(phi) result(r_hat)
  real(kind=real_kind), intent(in) :: phi(np,np)
  real(kind=real_kind)           :: r_hat(np,np)
  r_hat = elementwise_ij(phi, "r_hat_from_phi_scalar")
end function r_hat_from_phi_rank_ij
function r_hat_from_phi_rank_k(phi, k) result(r_hat)
  integer,              intent(in) :: k
  real(kind=real_kind), intent(in) :: phi(k)
  real(kind=real_kind)           :: r_hat(k)
  r_hat = elementwise_k(phi, k, "r_hat_from_phi_scalar")
end function r_hat_from_phi_rank_k
function r_hat_from_phi_no_rank(phi) result(r_hat)
  real(kind=real_kind), intent(in) :: phi
  real(kind=real_kind)           ::   r_hat
  r_hat =  r_hat_from_phi_scalar(phi)
end function r_hat_from_phi_no_rank
function z_from_phi_rank_ijk(phi, k) result(z)
  integer,              intent(in) :: k
  real(kind=real_kind), intent(in) :: phi(np,np,k)
  real(kind=real_kind)           :: z(np,np,k)
  z = elementwise_ijk(phi, k, "z_from_phi_scalar")
end function z_from_phi_rank_ijk
function z_from_phi_rank_ij(phi) result(z)
  real(kind=real_kind), intent(in) :: phi(np,np)
  real(kind=real_kind)           :: z(np,np)
  z = elementwise_ij(phi, "z_from_phi_scalar")
end function z_from_phi_rank_ij
function z_from_phi_rank_k(phi, k) result(z)
  integer,              intent(in) :: k
  real(kind=real_kind), intent(in) :: phi(k)
  real(kind=real_kind)           :: z(k)
  z = elementwise_k(phi, k, "z_from_phi_scalar")
end function z_from_phi_rank_k
function z_from_phi_no_rank(phi) result(z)
  real(kind=real_kind), intent(in) :: phi
  real(kind=real_kind)           :: z
  z =  z_from_phi_scalar(phi)
end function z_from_phi_no_rank
function g_from_z_rank_ijk(z, k) result(g)
  integer,              intent(in) :: k
  real(kind=real_kind), intent(in) :: z(np,np,k)
  real(kind=real_kind)           :: g(np,np,k)
  g = elementwise_ijk(z, k, "g_from_z_scalar")
end function g_from_z_rank_ijk
function g_from_z_rank_ij(z) result(g)
  real(kind=real_kind), intent(in) :: z(np,np)
  real(kind=real_kind)           :: g(np,np)
  g = elementwise_ij(z, "g_from_z_scalar")
end function g_from_z_rank_ij
function g_from_z_rank_k(z, k) result(g)
  integer,              intent(in) :: k
  real(kind=real_kind), intent(in) :: z(k)
  real(kind=real_kind)           :: g(k)
  g = elementwise_k(z, k, "g_from_z_scalar")
end function g_from_z_rank_k
function g_from_z_no_rank(z) result(g)
  real(kind=real_kind), intent(in) :: z
  real(kind=real_kind)           :: g
  g =  g_from_z_scalar(z)
end function g_from_z_no_rank
function g_from_phi_rank_ijk(phi, k) result(g)
  integer,              intent(in) :: k
  real(kind=real_kind), intent(in) :: phi(np,np,k)
  real(kind=real_kind)           :: g(np,np,k)
  g = elementwise_ijk(phi, k, "g_from_phi_scalar")
end function g_from_phi_rank_ijk
function g_from_phi_rank_ij(phi) result(g)
  real(kind=real_kind), intent(in) :: phi(np,np)
  real(kind=real_kind)           :: g(np,np)
  g = elementwise_ij(phi, "g_from_phi_scalar")
end function g_from_phi_rank_ij
function g_from_phi_rank_k(phi, k) result(g)
  integer,              intent(in) :: k
  real(kind=real_kind), intent(in) :: phi(k)
  real(kind=real_kind)           :: g(k)
  g = elementwise_k(phi, k, "g_from_phi_scalar")
end function g_from_phi_rank_k
function g_from_phi_no_rank(phi) result(g)
  real(kind=real_kind), intent(in) :: phi
  real(kind=real_kind)           :: g
  g =  g_from_phi_scalar(phi)
end function g_from_phi_no_rank
function dispatch(inval, funcname) result(outval)
  real(kind=real_kind), intent(in) :: inval
  character(*), intent(in) :: funcname
  real(kind=real_kind)          :: outval

  select case (funcname)
   case ("r_hat_from_phi_scalar")
     outval = r_hat_from_phi_scalar(inval)
   case ("z_from_phi_scalar")
     outval = z_from_phi_scalar(inval)
   case ("g_from_phi_scalar")
     outval = g_from_phi_scalar(inval)
   case ("g_from_z_scalar")
     outval = g_from_z_scalar(inval)
   end select
end function dispatch

function elementwise_ijk(inval, k, funcname) result(outval)
  integer,              intent(in) :: k
  real(kind=real_kind), intent(in) :: inval(np,np,k)
  real(kind=real_kind) :: outval(np,np,k)
  character (*), intent(in) :: funcname
  
  integer :: i, j, kk

  do kk=1,size(inval, 3)
    do j=1,size(inval, 2)
      do i=1,size(inval,1)
        outval(i,j,kk) = dispatch(inval(i, j, kk), funcname)
      end do
    end do
  end do
end function elementwise_ijk
function elementwise_ij(inval, funcname) result(outval)
  real(kind=real_kind), intent(in) :: inval(np,np)
  real(kind=real_kind) :: outval(np,np)
  character (*), intent(in) :: funcname  
  integer :: i, j

  do j=1,size(inval, 2)
    do i=1,size(inval,1)
      outval(i,j) = dispatch(inval(i, j), funcname)
    end do
  end do
end function elementwise_ij
function elementwise_k(inval, k, funcname) result(outval)
  integer,              intent(in) :: k
  real(kind=real_kind), intent(in) :: inval(k)
  real(kind=real_kind) :: outval(k)
  character (*), intent(in) :: funcname
  
  integer :: kk

  do kk=1,size(inval, 1)
    outval(kk) = dispatch(inval(kk), funcname)
  end do
end function elementwise_k


function g_from_z_scalar(z) result(g)
  real(kind=real_kind), intent(in) :: z
  real(kind=real_kind)             :: g
  if (atm_is_deep .and. .not. gravity_is_const) then
    g = gravit / ((rearth + z)/rearth)**2
  else
    g = gravit
  end if
end function g_from_z_scalar

function r_hat_from_phi_scalar(phi) result(r_hat)
  real(kind=real_kind), intent(in) :: phi
  real(kind=real_kind)             :: r_hat
  if (atm_is_deep) then 
    r_hat = (rearth + z_from_phi_scalar(phi))/rearth
  else
    r_hat = 1.0_real_kind
  end if

end function r_hat_from_phi_scalar

function z_from_phi_scalar(phi) result(z)
  real(kind=real_kind), intent(in) :: phi
  real(kind=real_kind)             :: z
  real(kind=real_kind)             :: b
  if (atm_is_deep .and. .not. gravity_is_const) then
    b = (2 * phi * rearth -gravit*rearth**2)
    z = -2*phi*rearth**2/(b-sqrt(b**2 - 4 * phi**2 * rearth**2) )
  else
    z = phi/gravit
  end if
end function z_from_phi_scalar


function g_from_phi_scalar(phi) result(g)
  real(kind=real_kind), intent(in) :: phi
  real(kind=real_kind)             :: g
  g = g_from_z_scalar(z_from_phi_scalar(phi))
end function g_from_phi_scalar


function quasi_hydrostatic_terms(u_top, lat, phi_top) result(terms)
  real(kind=real_kind), intent(in) :: u_top(np,np,2)
  real(kind=real_kind), intent(in) :: lat(np,np)
  real(kind=real_kind), intent(in) :: phi_top(np,np)
  real(kind=real_kind) :: terms(np,np)
  terms = - (u_top(:,:,1)**2 + u_top(:,:,2)**2)/(rearth + z_from_phi(phi_top)) - 2 * omega * cos(lat) * u_top(:,:,1) 

end function quasi_hydrostatic_terms
end module 
