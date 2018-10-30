  subroutine build_grad_mass(c,e,z,nr,r,vr,ia,ie)
! Calculated the gradient of the scalar relativistic mass
! grad_k ln(M(r)) = - 2 M(r)*v_soc(r)*r*(e_r)_k
! with
! vsoc = 1/(2Mc)^2 1/r dV/dr
! M = m + (E - V(r))/(2c^2)
! m = 1/2, e^2 = 2 in Ry units
! e-e potential     is repulsive  -> positive
! nuclear potential is attractive -> negative
  use global, only: lmmax, lm2i, i4b, r8b, c8b, grad_mass, nasusc, npanat, ircutat
  use mod_derivative_panels

  implicit none

! Speed of light
  real(kind=r8b),    intent(in)  :: c
! Complex energy
  complex(kind=c8b), intent(in)  :: e
! Atomic number
  real(kind=r8b),    intent(in)  :: z
! Energy point number
  integer(kind=i4b), intent(in)  :: ie
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
  complex(kind=c8b) :: mass(nr), work(nr)
  real(kind=r8b)    :: dvdr(nr)
  integer(kind=i4b) :: ir,i
  integer(kind=i4b) :: ilmxyz(1:3)
! ilm for unit vector in direction x,y,z
  ilmxyz(1)=lm2i(1,1);ilmxyz(2)=lm2i(-1,1);ilmxyz(3)=lm2i(0,1)
! Relativistic mass (without c, differs from build_vsoc)
  mass(1:nr) = 1.d0 + (e - vr(1:nr) + 2.d0*z/r(1:nr))/c**2
! write out mass
  write(229,'("# r mass           for ia,ie=",2i4)') ia,ie
  do ir=1,nr
    write(229,'(10e18.9)') r(ir), mass(ir)
  end do
!! dV/dr for e-e potential
!  call calc_derivative_panels(vr(1:nr),dvdr(1:nr),r,nr,npanat(ia),ircutat(:,ia))
!! 1/r dV/dr
!  vsoc(1:nr) = 2.d0*z/(r(1:nr)**3) + dvdr/r(1:nr)
!! multiply by the inverse relativistic mass and the speed of light
!  vsoc(1:nr) = vsoc(1:nr)/(mass(1:nr))**2/c**2
!! calculate grad_mass=-mass(r)*vsoc*r
!  do i=1,3
!    grad_mass(i,1:nr,ilmxyz(i),ia,ie)=-1.d0*mass(1:nr)*vsoc(1:nr)*r(1:nr)
!  end do
! calculate d lnM/ dr directly by finite differences  
  call calc_derivative_panels(log(mass(1:nr)),work(1:nr),r,nr,npanat(ia),ircutat(:,ia))
  do i=1,3
    grad_mass(i,1:nr,ilmxyz(i),ia,ie) = work(1:nr)
  end do
  write(239,'("# r grad_mass           for ia,ie=",2i4)') ia,ie
  do ir = 1,nr
    write(239,'(100e18.9)') r(ir), (grad_mass(i,ir,ilmxyz(i),ia,ie),i=1,3)
  end do
! All done!
  end subroutine build_grad_mass
