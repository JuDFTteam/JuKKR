!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_corehff
  
  private
  public :: corehff

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculate rel. hyberfine fields of core state
  !> Author: 
  !> Category: KKRhost, core-electrons, physical-observables
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> CALCULATE THE RELATIVISTIC HYPERFINEFIELDS FOR THE
  !>                CURRENT  CORE STATE S              
  !>
  !> THE WAVE FUNCTION  {G(K,S),F(K,S)}  IS NORMALIZED TO 1
  !-------------------------------------------------------------------------------
  subroutine corehff(kap1, kap2, mj, s, nsol, bhf, gck, fck, rc, drdic, rnuc, nzero, nrc)

    use :: mod_ylag, only: ylag
    use :: mod_rint4pts, only: rint4pts
    use :: mod_datatypes, only: dp
    implicit none

    ! PARAMETER definitions
    real (kind=dp), parameter :: eps = 1.0e-12_dp
    real (kind=dp) :: e0, a0, cautog

    ! CONVERSION FACTOR FOR HYPERFINE FIELDS FROM A.U. TO GAUSS
    ! ELECTRON CHARGE     IN ESU
    ! BOHR-RADIUS         IN CM

    parameter (e0=1.6021892e-19_dp*2.997930e+09_dp, a0=0.52917706e-08_dp, cautog=e0/(a0*a0))

    ! Dummy arguments
    integer :: kap1, kap2, nrc, nsol, nzero, s
    real (kind=dp) :: mj, rnuc
    real (kind=dp) :: bhf(2, 2), drdic(nrc), fck(2, 2, nrc), gck(2, 2, nrc), rc(nrc)

    ! Local variables
    real (kind=dp) :: ame(2, 2), xx(5), yi(nrc), yy(5), zi(nrc)
    integer :: i, k1, k2, n

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! ANGULAR HYPERFINE MATRIX ELEMENTS   SEE E.G.  E.M.ROSE
    ! THE FACTOR  I  HAS BEEN OMITTED
    ame(1, 1) = 4.0e0_dp*kap1*mj/(4.0e0_dp*kap1*kap1-1.0e0_dp)
    if (nsol==2) then
      ame(2, 2) = 4.0e0_dp*kap2*mj/(4.0e0_dp*kap2*kap2-1.0e0_dp)
      ame(1, 2) = sqrt(0.25e0_dp-(mj/real(kap1-kap2,kind=dp))**2)
      ame(2, 1) = ame(1, 2)
    end if
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    do k1 = 1, nsol
      do k2 = 1, nsol
        do n = 1, nzero
          yi(n) = drdic(n)*(gck(k1,s,n)*fck(k2,s,n)+fck(k1,s,n)*gck(k2,s,n))
        end do
        call rint4pts(yi, nzero, zi)
        if (abs(rnuc)>eps) then
          do i = 1, 5
            xx(i) = rc(nzero-5+i)
            yy(i) = zi(nzero-5+i)
          end do
          zi(nzero) = ylag(rnuc, xx, yy, 0, 4, 5)
        end if
        xx(1) = 1.0e0_dp
        ! !RC( 1)
        xx(2) = 6.0e0_dp
        ! !RC( 6)
        xx(3) = 11.0e0_dp
        ! !RC(11)
        yy(1) = zi(nzero) - zi(1)
        yy(2) = zi(nzero) - zi(6)
        yy(3) = zi(nzero) - zi(11)
        bhf(k1, k2) = cautog*ame(k1, k2)*ylag(0.0e0_dp, xx, yy, 0, 2, 3)
      end do
    end do

  end subroutine corehff

end module mod_corehff
