!-------------------------------------------------------------------------------
! SUBROUTINE: RHOVAL0
!> @note Jonathan Chico Apr. 2019: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
subroutine rhoval0(ez, drdi, rmesh, ipan, ircut, irws, thetas, dos0, dos1, &
  irm, lmax)
!
  use :: constants
  use :: global_variables

  implicit none

! .. Input variables
  integer, intent (in) :: irm !< Maximum number of radial points
  integer, intent (in) :: lmax !< Maximum l component in wave function expansion
  integer, intent (in) :: ipan !< Number of panels in non-MT-region
  integer, intent (in) :: irws !< R point at WS radius
  double complex, intent (in) :: ez
  integer, dimension (0:ipand), intent (in) :: ircut !< R points of panel borders
  double precision, dimension (irm), intent (in) :: drdi !< Derivative dr/di
  double precision, dimension (irm), intent (in) :: rmesh
  double precision, dimension (irid, nfund), intent (in) :: thetas !< shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
! .. Output variables
  double complex, intent (out) :: dos0
  double complex, intent (out) :: dos1
! .. Local Scalars
  integer :: ir, l, l1, imt1
  integer :: lmaxd1
  double precision :: c0ll
  double complex :: ek, ciek, denl
! .. Local Arrays ..
  double complex, dimension (0:lmax+1) :: bessjw
  double complex, dimension (0:lmax+1) :: bessyw
  double complex, dimension (0:lmax+1) :: hankws
  double complex, dimension (irm, 0:lmax) :: pz
  double complex, dimension (irm, 0:lmax) :: qz
  double complex, dimension (irm, 0:lmax+1) :: cden0
  double complex, dimension (irm, 0:lmax+1) :: cden1

! .. Intrinsic Functions
  intrinsic :: atan, dble, sqrt
!
  lmaxd1 = lmax + 1
  ek = sqrt(ez)
  c0ll = 1.0d0/sqrt(16.0d0*atan(1.0d0))
  ciek = ci*ek
!
!----------------------------------------------------------------------------
  do ir = 2, irws
    call beshan(hankws, bessjw, bessyw, rmesh(ir)*ek, lmaxd1)
    do l = 0, lmax
      pz(ir, l) = bessjw(l)*rmesh(ir)
      qz(ir, l) = (bessyw(l)-ci*bessjw(l))*rmesh(ir)
    end do
  end do
  imt1 = ircut(1)
  do l1 = 0, lmaxd1
    cden0(1, l1) = czero
    cden1(1, l1) = czero
  end do
  do ir = 2, irws
    cden0(ir, 0) = ek*pz(ir, 0)*qz(ir, 0)
    cden1(ir, 0) = ek*pz(ir, 0)**2*(0.d0, -1.d0)
    cden1(ir, lmaxd1) = ciek*rmesh(ir)**2
  end do
  do l1 = 1, lmax
    do ir = 2, irws
      cden0(ir, l1) = ek*pz(ir, l1)*qz(ir, l1)*(l1+l1+1)
      cden1(ir, l1) = ek*pz(ir, l1)**2*(0.d0, -1.d0)*(l1+l1+1)
    end do
  end do
!
  do l1 = 0, lmaxd1 !LMAXD1
    if (ipan>1) then
      do ir = imt1 + 1, irws
        cden0(ir, l1) = cden0(ir, l1)*thetas(ir-imt1, 1)*c0ll
        cden1(ir, l1) = cden1(ir, l1)*thetas(ir-imt1, 1)*c0ll
      end do
    end if
    call csimpk(cden0(1,l1), denl, ipan, ircut, drdi)
    dos0 = dos0 + denl
    call csimpk(cden1(1,l1), denl, ipan, ircut, drdi)
    dos1 = dos1 + denl
  end do
!
end subroutine
