  subroutine build_grad_mass(c,e,z,nr,r,vr,ia)
! Calculated the gradient of the scalar relativistic mass
! grad_k ln(M(r)) = - 2 M(r)*v_soc(r)*r*(e_r)_k
! with
! vsoc = 1/(2Mc)^2 1/r dV/dr
! M = m + (E - V(r))/(2c^2)
! m = 1/2, e^2 = 2 in Ry units
! e-e potential     is repulsive  -> positive
! nuclear potential is attractive -> negative
  use global, only: lmmax, lm2i, i4b, r8b, c8b, grad_mass, nasusc, npanat, ircutat
  use derivative_panels

  implicit none

! Speed of light
  real(kind=r8b),    intent(in)  :: c
! Complex energy
  complex(kind=c8b), intent(in)  :: e
! Atomic number
  real(kind=r8b),    intent(in)  :: z
! Atom number
  integer(kind=i4b), intent(in)  :: ia
! Number of radial points
  integer(kind=i4b), intent(in)  :: nr
! Radial mesh values
  real(kind=r8b),    intent(in)  :: r(nr)
! Spherical spin-averaged radial potential
  real(kind=r8b),    intent(in)  :: vr(nr)
! -------------------------------------------------------------------
! Complex SOC potential
  complex(kind=c8b) :: vsoc(nr)
  complex(kind=c8b) :: mass(nr)
  real(kind=r8b)    :: dvdr(nr)
  integer(kind=i4b) :: ir,i
  integer(kind=i4b) :: ilmxyz(1:3)
! ilm for unit vector in direction x,y,z
  ilmxyz(1)=lm2i(1,1);ilmxyz(2)=lm2i(-1,1);ilmxyz(3)=lm2i(0,1)
! Relativistic mass multiplied by c
  mass(1:nr) = c + (e - vr(1:nr) + 2.d0*z/r(1:nr))/c
! dV/dr for e-e potential
  call calc_derivative_panels(vr(1:nr),dvdr(1:nr),r,nr,npanat(ia),ircutat(:,ia))
!! forward difference
!  dvdr(1) = (vr(2) - vr(1))/(r(2) - r(1))
!! backward difference
!  dvdr(nr) = (vr(nr) - vr(nr-1))/(r(nr) - r(nr-1))
!! centered differences
!  do ir=2,nr-1
!    dvdr(ir) = 0.5d0*(vr(ir+1) - vr(ir))/(r(ir+1) - r(ir))
!    dvdr(ir) = dvdr(ir) + 0.5d0*(vr(ir) - vr(ir-1))/(r(ir) - r(ir-1))
!  end do
! 1/r dV/dr
  vsoc(1:nr) = 2.d0*z/(r(1:nr)**3) + dvdr/r(1:nr)
! multiply by the inverse relativistic mass
  vsoc(1:nr) = vsoc(1:nr)/(mass(1:nr))**2
! calculate grad_mass=-2*mass(r)*vsoc*r
  do i=1,3
    grad_mass(i,1:nr,ilmxyz(i),ia)=-2.d0*mass(1:nr)*vsoc(1:nr)*r(1:nr)
  end do
  write(*,'(" grad_mass constructed for atom ",i4)') ia
  write(*,'(" r grad_mass x/y/z")')
  do ir = 1,nr
    write(*,'(100e18.9)') r(ir), (grad_mass(i,ir,ilmxyz(i),ia),i=1,3)
  end do
! All done!
  end subroutine build_grad_mass
