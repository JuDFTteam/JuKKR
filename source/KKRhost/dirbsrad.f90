!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_dirbsrad

contains

  !-------------------------------------------------------------------------------
  !> Summary: Derivative of coupled radial Dirac equations 
  !> Author: 
  !> Category: KKRhost, dirac, single-site
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Supply the derivatives for the coupled set of
  !> radial dirac equation in case of a spin-dependent potential
  !-------------------------------------------------------------------------------
  !> @note Jonathan Chico: This makes use of common blocks, one should try to remove this
  !-------------------------------------------------------------------------------
  subroutine dirbsrad(xbs, y, dydx, drdi, b, v, r, nmesh)

    use :: mod_datatypes, only: dp
    use :: mod_dirbslag, only: dirbslag
    implicit none

    include 'sprkkr_rmesh.dim'


    ! PARAMETER definitions
    integer :: nlag
    parameter (nlag=3)

    ! COMMON variables
    real (kind=dp) :: cgd(2), cgmd(2), cgo, kap(2)
    real (kind=dp) :: csqr
    complex (kind=dp) :: ebs
    integer :: nradbs, nsolbs
    common /commbs/ebs, csqr, cgd, cgmd, cgo, kap, nsolbs, nradbs

    ! Dummy arguments
    integer :: nmesh
    real (kind=dp) :: xbs
    real (kind=dp) :: b(nrmax), drdi(nrmax), r(nrmax), v(nrmax)
    complex (kind=dp) :: dydx(ncfmax), y(ncfmax)

    ! Local variables
    real (kind=dp) :: bbs, bpp, bqq, drdibs, rbs, vbs
    complex (kind=dp) :: emvpp, emvqq
    integer :: i, j, k, m

    call dirbslag(xbs, vbs, bbs, rbs, drdibs, v, b, r, drdi, nradbs, nlag, nmesh)

    emvqq = (ebs-vbs+csqr)*drdibs/csqr
    emvpp = -(ebs-vbs)*drdibs
    bqq = bbs*drdibs/csqr
    bpp = bbs*drdibs


    m = 0
    do j = 1, nsolbs
      do i = 1, nsolbs
        k = 3 - 4*(i-1)
        dydx(m+1) = -kap(i)*y(m+1)/rbs*drdibs + (emvqq+bqq*cgmd(i))*y(m+2)
        dydx(m+2) = kap(i)*y(m+2)/rbs*drdibs + (emvpp+bpp*cgd(i))*y(m+1) + bpp*cgo*y(m+k)
        m = m + 2
      end do
    end do
  end subroutine dirbsrad

end module mod_dirbsrad
