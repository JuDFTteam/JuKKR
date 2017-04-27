  subroutine build_vsoc(c,e,z,nr,r,vr,socscaling,vsoc,numpan,numrcut)
! Assembles the radial part of the SOC potential
! vsoc = 1/(2Mc)^2 1/r dV/dr
! M = m + (E - V(r))/(2c^2)
! m = 1/2, e^2 = 2 in Ry units
! e-e potential     is repulsive  -> positive
! nuclear potential is attractive -> negative
! Modified for npan > 1
! Same radial shift for all atoms 
  use global, only: i4b, r8b, c8b, rmesh, nrpts0

  implicit none

! Speed of light
  real(kind=r8b),    intent(in)  :: c
! Complex energy
  complex(kind=c8b), intent(in)  :: e
! Atomic number
  real(kind=r8b),    intent(in)  :: z
! Number of radial points
  integer(kind=i4b), intent(in)  :: nr
! Radial mesh values
  real(kind=r8b),    intent(in)  :: r(nr)
! Spherical spin-averaged radial potential
  real(kind=r8b),    intent(in)  :: vr(nr)
! Scaling of potential strength
  real(kind=r8b),    intent(in)  :: socscaling
! Complex SOC potential
  complex(kind=c8b), intent(out) :: vsoc(nr)
! --> Number of panels > 1
  integer(kind=i4b), intent(in) :: numpan, numrcut(numpan+1) 
! -------------------------------------------------------------------
  complex(kind=c8b) :: mass(nr)
  real(kind=r8b)    :: dvdr(nr)
  integer(kind=i4b) :: ir, ip, ist, ien

! Relativistic mass multiplied by c
  mass(1:nr) = c + (e - vr(1:nr) + 2.d0*z/r(1:nr))/c

  dvdr = 0.d0 
 
! Loop over panels
  do ip = 1, numpan

    ! Begin and end of panels
    if (ip == 1) then
      ist = numrcut(ip) + 1            ! Start from 1  
    else
      ist = numrcut(ip) - nrpts0(1) + 2                                
    endif
      ien = numrcut(ip+1) - nrpts0(1) + 1                                

    ! dV/dr for e-e potential
    ! forward difference 
    ! Carful the first point of the panel is repeated
    dvdr(ist) = (vr(ist+1) - vr(ist))/(r(ist+1) - r(ist))
 
    ! backward differences
    dvdr(ien) = (vr(ien) - vr(ien-1))/(r(ien) - r(ien-1))

    ! centered differences
    do ir =  ist+1, ien-1 
      dvdr(ir) = 0.5d0*(vr(ir+1) - vr(ir))/(r(ir+1) - r(ir))
      dvdr(ir) = dvdr(ir) + 0.5d0*(vr(ir) - vr(ir-1))/(r(ir) - r(ir-1))
    end do

  end do ! panels
 
! 1/r dV/dr
  vsoc(1:nr) = 2.d0*z/(r(1:nr)**3) + dvdr/r(1:nr)
! multiply by the inverse relativistic mass
  vsoc(1:nr) = socscaling*vsoc(1:nr)/(mass(1:nr))**2
 
! All done!
  end subroutine build_vsoc
