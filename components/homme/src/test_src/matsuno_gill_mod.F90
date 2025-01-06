#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module matsuno_gill_mod


  use coordinate_systems_mod, only: spherical_polar_t
  use dimensions_mod,         only: nlev,np,qsize,nlevp
  use element_mod,            only: element_t
  use element_state,          only: timelevels
  use element_ops,            only: set_thermostate, get_temperature
  use hybrid_mod,             only: hybrid_t
  use hybvcoord_mod,          only: hvcoord_t
  use kinds,                  only: real_kind, iulog
  use physical_constants,     only: p0, kappa,g, dd_pi, Rgas,Cp
  use physics_mod,            only: prim_condense
  use time_mod,               only: secpday
#ifndef HOMME_WITHOUT_PIOLIBRARY
  use common_io_mod,          only: infilenames
#endif

implicit none
private

  real (kind=real_kind), public, parameter :: sigma_b = 0.70D0
  real (kind=real_kind), public, parameter :: k_Tv     = 1.0D0/(40.0D0*secpday)
  real (kind=real_kind), public, parameter :: bv_sq   = 1.0D-4
  real (kind=real_kind), public, parameter :: T0 = 300
  real (kind=real_kind), public, parameter :: q0   = 5d0 /secpday
  real (kind=real_kind), public, parameter :: dlam   = 1.0D0 * 30D0 * dd_pi/180D0 
  real (kind=real_kind), public, parameter :: dphi   = 1.0D0 * 10D0 * dd_pi/180D0 

  public :: mg_v_forcing
  public :: mg_T_forcing
  public :: mg_init_state
  public :: mg_forcing

contains

  subroutine mg_forcing(elemin,hvcoord,nm1,nm1_Q,dt)

    type (element_t)                 :: elemin
    type (hvcoord_t)                 :: hvcoord
    integer                          :: nm1,nm1_Q  ! timelevel to use
    real (kind=real_kind)                   :: dt

    ! local
    real (kind=real_kind)                   :: pmid,r0,r1,dtf_q,dp,rdp,FQ
    real (kind=real_kind), dimension(np,np) :: psfrc 
    real (kind=real_kind)                   :: temperature(np,np,nlev)
    real (kind=real_kind)                   :: v(np,np,3,nlev)
    real (kind=real_kind)                   :: fv(np,np,3,nlev)
    integer                                 :: i,j,k,q

    dtf_q = dt
    call get_temperature(elemin,temperature,hvcoord,nm1)
        
    do j=1,np
       do i=1,np
          psfrc(i,j) = (elemin%state%ps_v(i,j,nm1))
       end do
    end do

    elemin%derived%FT(:,:,:) = elemin%derived%FT(:,:,:) + &
         mg_T_forcing(hvcoord,psfrc(1,1),               &
         temperature,elemin%spherep,np, nlev)

    v(:,:,1:2,:) = elemin%state%v(:,:,1:2,:,nm1)
#if ( defined MODEL_THETA_L ) 
    v(:,:,3,:) = elemin%state%w_i(:,:,1:nlev,nm1)  ! dont apply at surface
#else
    v(:,:,3,:) = 0
#endif

    fv = mg_v_forcing(hvcoord,psfrc(1,1),v,np,nlev)

#if ( defined MODEL_THETA_L ) 
    elemin%derived%FM(:,:,1:3,:) = elemin%derived%FM(:,:,1:3,:) + fv(:,:,1:3,:)
#else
    elemin%derived%FM(:,:,1:2,:) = elemin%derived%FM(:,:,1:2,:) + fv(:,:,1:2,:)
#endif

    if (qsize>=1) then
       ! HS with tracer  (Galewsky type forcing, with flux of  2.3e-5 kg/m^2/s
       ! MASS in kg/m^2   = < Q dp_in_Pa / g >   
       ! flux in kg/m^2/s = < FQ dp_in_Pa / g >   
       ! We want < FQ dp_in_Pa / g > = 2.3e-5  so:  FQ = 2.3e-5*g/dp_in_Pa 

       ! lowest layer thickness, in Pa
       dp = ( hvcoord%hyai(nlev+1) - hvcoord%hyai(nlev) ) + &
               ( hvcoord%hybi(nlev+1) - hvcoord%hybi(nlev) )*1000*100
       rdp = 1./ dp
       q=1
       do j=1,np
          do i=1,np
             FQ = rdp * g * 2.3E-5 * COS(elemin%spherep(i,j)%lat)**2
             elemin%derived%FQ(i,j,nlev,q) =elemin%derived%FQ(i,j,nlev,q)+FQ
          enddo
       enddo

       do j=1,np
          do i=1,np
             do k=1,nlev
                pmid = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*(elemin%state%ps_v(i,j,nm1))
                r0=elemin%state%Q(i,j,k,q)
                r1=r0
                call Prim_Condense(r1,temperature(i,j,k),pmid)
                elemin%derived%FQ(i,j,k,q) = elemin%derived%FQ(i,j,k,q) + &
                     (r1-r0)/(dtf_q)
             enddo
          enddo
       enddo
    endif

  end subroutine mg_forcing
  function z_given_p(p) result(z)
      real (kind=real_kind), intent(in) :: p
      real (kind=real_kind):: z
      z = (-g/bv_sq)*log( (T0/(g**2/(bv_sq*Cp)))*( (p/p0)**(Rgas/Cp) - 1.0  ) + 1.0 )
      !z = (-g/Nsq)*np.log( (T0/bigG)*( (p/p0)**(Rd/Cp) - 1.0  ) + 1.0 )
  end function 
  function T_given_p(p) result(T)
      real (kind=real_kind), intent(in) :: p
      real (kind=real_kind) :: T
      T = T0 * (p/p0)**kappa * exp(bv_sq/g * z_given_p(p))
  end function
 
  function mg_v_forcing(hvcoord,ps,v,npts,nlevels) result(mg_v_frc)

    integer, intent(in)               :: npts
    integer, intent(in)               :: nlevels
    type (hvcoord_t), intent(in)       :: hvcoord
    real (kind=real_kind), intent(in) :: ps(npts,npts)

    real (kind=real_kind), intent(in) :: v(npts,npts,3,nlevels)
    real (kind=real_kind)             :: mg_v_frc(npts,npts,3,nlevels)

    ! Local variables

    integer i,j,k
    real (kind=real_kind) :: k_v
    real (kind=real_kind) :: p,etam

    do k=1,nlevels
       do j=1,npts
          do i=1,npts
             p    = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*ps(i,j)
             mg_v_frc(i,j,1,k) = -k_Tv*v(i,j,1,k)
             mg_v_frc(i,j,2,k) = -k_Tv*v(i,j,2,k)

             mg_v_frc(i,j,3,k) = -k_Tv*v(i,j,3,k)
          end do
       end do
    end do

  end function mg_v_forcing

  function mg_T_forcing(hvcoord,ps,T,sphere,npts,nlevels) result(mg_T_frc)

    integer, intent(in) :: npts
    integer, intent(in) :: nlevels

    type (hvcoord_t), intent(in)          :: hvcoord
    real (kind=real_kind), intent(in)    :: ps(npts,npts)
    real (kind=real_kind), intent(in)    :: T(npts,npts,nlevels)
    type (spherical_polar_t), intent(in) :: sphere(npts,npts)

    real (kind=real_kind)             :: mg_T_frc(npts,npts,nlevels)

    ! Local variables

    real (kind=real_kind) :: p,logprat,pratk,Teq,T_forcing
    real (kind=real_kind) :: logps0,etam
    real (kind=real_kind) :: lat,snlat

    real (kind=real_kind) :: k_t(npts,npts)
    real (kind=real_kind) :: snlatsq(npts,npts)
    real (kind=real_kind) :: latlon_comp(npts,npts), z_top(npts,npts)
    real (kind=real_kind) :: cslatsq(npts,npts)

    real (kind=real_kind) :: rec_one_minus_sigma_b

    integer i,j,k

    logps0    = LOG(hvcoord%ps0 )

    do j=1,npts
       do i=1,npts
         snlat        = SIN(sphere(i,j)%lat)
         if (abs(sphere(i,j)%lon-dd_pi/2) < dlam .and. abs(sphere(i,j)%lat) < dphi ) then
         latlon_comp(i,j) = cos((sphere(i,j)%lon-dd_pi/2)/(2 * dlam))**2 * cos(sphere(i,j)%lat/(2*dphi))**2
         else
           latlon_comp(i,j) = 0
         end if
         z_top(i,j) = z_given_p(hvcoord%hyai(1)*hvcoord%ps0 + hvcoord%hybi(1)*ps(i,j))
       end do
    end do

    rec_one_minus_sigma_b = 1.0D0/(1.0D0 - sigma_b)

    do k=1,nlevels
       do j=1,npts
          do i=1,npts
             p         = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*ps(i,j)

             T_forcing = q0 * latlon_comp(i,j) * sin(dd_pi * z_given_p(p)/z_top(i,j))
#if 0
             ! ======================================
             ! This is a smooth forcing 
             ! for debugging purposes only...
             ! ======================================

#endif
             mg_T_frc(i,j,k)= -k_Tv*(T(i,j,k)-T_given_p(p)) + T_forcing
          end do
       end do
    end do
     
  end function mg_T_forcing

  subroutine mg_init_state(elem, hybrid, hvcoord,nets,nete)

    type(element_t),        intent(inout) :: elem(:)
    type(hybrid_t),         intent(in)    :: hybrid                   ! hybrid parallel structure
    type (hvcoord_t),       intent(in)    :: hvcoord
    integer,                intent(in)    :: nets
    integer,                intent(in)    :: nete

    ! Local variables
    
    integer ie,i,j,k,q,tl
    integer :: nm1 
    integer :: n0 
    integer :: np1
    real (kind=real_kind) :: lat_mtn,lon_mtn,r_mtn,h_mtn,rsq,lat,lon,p_tmp
    real (kind=real_kind) :: temperature(np,np,nlev),p(np,np),exner(np,np),ps(np,np)

    if (hybrid%masterthread) write(iulog,*) 'initializing Held-Suarez primitive equations test'

    nm1= 1
    n0 = 2
    np1= 3

    do ie=nets,nete

       elem(ie)%state%ps_v(:,:,n0) =p0
       elem(ie)%state%ps_v(:,:,nm1)=p0
       elem(ie)%state%ps_v(:,:,np1)=p0

       elem(ie)%state%v(:,:,:,:,n0) =0.0D0
       elem(ie)%state%v(:,:,:,:,nm1)=elem(ie)%state%v(:,:,:,:,n0)
       elem(ie)%state%v(:,:,:,:,np1)=elem(ie)%state%v(:,:,:,:,n0)

#ifdef MODEL_THETA_L
       elem(ie)%state%w_i = 0.0
#endif


       ! if topo file was given in the namelist, PHIS was initilized in prim_main
       ! otherwise assume 0
#ifndef HOMME_WITHOUT_PIOLIBRARY
       if (infilenames(1)=='') then
          elem(ie)%state%phis(:,:)=0.0D0
       endif
#endif

       do k=1,nlev
         p_tmp = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem(ie)%state%ps_v(1,1,n0)
         temperature(1, 1, k) = T_given_p(p_tmp)
       end do
       do i=1,np
         do j=1,np
           temperature(i,j, :) = temperature(1,1, :)
         end do
       end do




       !if (qsize>=1) then
       !q=1
       !elem(ie)%state%Q(:,:,:,q) =0  ! moist HS tracer IC=0
       !do q=2,qsize
       !   elem(ie)%state%Q(:,:,:,q) =temperature(:,:,:)/400
       !enddo
       !endif
       ps=elem(ie)%state%ps_v(:,:,n0)
       call set_thermostate(elem(ie),ps,temperature,hvcoord) ! sets z based on hydrostatic presumption
    end do


  end subroutine mg_init_state


end module matsuno_gill_mod

