!------------------------------------------------------------------------------------
!> Summary: Calculate the contribution to the DOS from free space
!> Author: 
!> Calculate the contribution to the DOS from free space
!------------------------------------------------------------------------------------
module mod_rhoval0

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculate the contribution to the DOS from free space 
  !> Author: 
  !> Category: physical-observables, KKRhost
  !> Deprecated: False 
  !> Calculate the contribution to the DOS from free space
  !-------------------------------------------------------------------------------
  subroutine rhoval0(ez,drdi,rmesh,ipan,ircut,irws,thetas,dos0,dos1,irm,lmax)

    use :: mod_constants
    use :: global_variables
    use :: mod_datatypes, only: dp
    use :: mod_beshan
    use :: mod_csimpk

    implicit none

    ! .. Input variables
    integer, intent (in) :: irm    !! Maximum number of radial points
    integer, intent (in) :: lmax   !! Maximum l component in wave function expansion
    integer, intent (in) :: ipan   !! Number of panels in non-MT-region
    integer, intent (in) :: irws   !! R point at WS radius
    complex (kind=dp), intent (in) :: ez
    integer, dimension (0:ipand), intent (in) :: ircut !! R points of panel borders
    real (kind=dp), dimension (irm), intent (in) :: drdi !! Derivative dr/di
    real (kind=dp), dimension (irm), intent (in) :: rmesh
    real (kind=dp), dimension (irid, nfund), intent (in) :: thetas !! shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
    ! .. Output variables
    complex (kind=dp), intent (out) :: dos0
    complex (kind=dp), intent (out) :: dos1
    ! .. Local Scalars
    integer :: ir, l, l1, imt1
    integer :: lmaxd1
    real (kind=dp) :: c0ll
    complex (kind=dp) :: ek, ciek, denl
    ! .. Local Arrays ..
    complex (kind=dp), dimension (0:lmax+1) :: bessjw
    complex (kind=dp), dimension (0:lmax+1) :: bessyw
    complex (kind=dp), dimension (0:lmax+1) :: hankws
    complex (kind=dp), dimension (irm, 0:lmax) :: pz
    complex (kind=dp), dimension (irm, 0:lmax) :: qz
    complex (kind=dp), dimension (irm, 0:lmax+1) :: cden0
    complex (kind=dp), dimension (irm, 0:lmax+1) :: cden1

    ! .. Intrinsic Functions
    intrinsic :: atan, sqrt

    lmaxd1 = lmax + 1
    ek = sqrt(ez)
    c0ll = 1.0e0_dp/sqrt(16.0e0_dp*atan(1.0e0_dp))
    ciek = ci*ek

    ! initialize to zero
    cden0(:, :) = czero
    cden1(:, :) = czero
    dos0 = czero
    dos1 = czero

    ! ----------------------------------------------------------------------------
    do ir = 2, irws
      call beshan(hankws, bessjw, bessyw, rmesh(ir)*ek, lmaxd1)
      do l = 0, lmax
        pz(ir, l) = bessjw(l)*rmesh(ir)
        qz(ir, l) = (bessyw(l)-ci*bessjw(l))*rmesh(ir)
      end do
    end do
    imt1 = ircut(1)
    do ir = 2, irws
      cden0(ir, 0) = ek*pz(ir, 0)*qz(ir, 0)
      cden1(ir, 0) = ek*pz(ir, 0)**2*(0.e0_dp, -1.e0_dp)
      cden1(ir, lmaxd1) = ciek*rmesh(ir)**2
    end do
    do l1 = 1, lmax
      do ir = 2, irws
        cden0(ir, l1) = ek*pz(ir, l1)*qz(ir, l1)*(l1+l1+1)
        cden1(ir, l1) = ek*pz(ir, l1)**2*(0.e0_dp, -1.e0_dp)*(l1+l1+1)
      end do
    end do

    do l1 = 0, lmaxd1              ! LMAXD1
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

  end subroutine rhoval0

end module mod_rhoval0
