!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Transformation of the wavefunctions for non spherical potentials.
!> Author: B. Drittler
!> To determine the non-spherical wavefunctions the potential has to be lm1 and lm2 dependent.
!> the potential is stored only as lm dependent , therefore a transformation in the
!> following way has to be done:
!> \begin{equation}
!> vnsll(r,lm1,lm2) = \sum_{lm3} \left\{  c(lm1,lm2,lm3) vins(r,lm3)\right\}
!> \end{equation}
!> where \(c(lm1,lm2,lm3)\) are the gaunt coeffients. (see notes by B. Drittler)
!------------------------------------------------------------------------------------
!> @note attention : The gaunt coeffients are stored in an index array only for lm1.gt.lm2
!> (see subroutine gaunt)
!> R. Zeller Sep. 2000: modified
!> @endnote
!------------------------------------------------------------------------------------
module mod_vllns
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Transformation of the wavefunctions for non spherical potentials.
  !> Author: B. Drittler
  !> Category: special-functions, potential, KKRhost 
  !> Deprecated: False
  !> To determine the non-spherical wavefunctions the potential has to be lm1 and lm2 dependent.
  !> the potential is stored only as lm dependent , therefore a transformation in the
  !> following way has to be done :
  !> \begin{equation}
  !> vnsll(r,lm1,lm2) = \sum_{lm3} \left\{  c(lm1,lm2,lm3) vins(r,lm3)\right\}
  !> \end{equation}
  !> where \(c(lm1,lm2,lm3)\) are the gaunt coeffients. (see notes by B. Drittler)
  !-------------------------------------------------------------------------------
  !> @note attention : The gaunt coeffients are stored in an index array only for lm1.gt.lm2
  !> (see subroutine gaunt)
  !> R. Zeller Sep. 2000: modified
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine vllns(vnspll,vins,cleb,icleb,iend,irm,ncleb,lmpot,irmind,lmmaxd)

    implicit none

    ! .. Input variables
    integer, intent (in) :: irm    !! Maximum number of radial points
    integer, intent (in) :: iend
    integer, intent (in) :: ncleb  !! Number of Clebsch-Gordon coefficients
    integer, intent (in) :: lmpot  !! (LPOT+1)**2
    integer, intent (in) :: irmind !! IRM-IRNSD
    integer, intent (in) :: lmmaxd !! (KREL+KORBIT+1)(LMAX+1)^2
    ! .. Array Arguments
    integer, dimension (ncleb, 4), intent (in) :: icleb !! Pointer array
    real (kind=dp), dimension (ncleb, 2), intent (in) :: cleb !! GAUNT coefficients (GAUNT)
    real (kind=dp), dimension (irmind:irm, lmpot), intent (in) :: vins !! Non-spherical part of the potential
    ! .. Output variables
    real (kind=dp), dimension (lmmaxd, lmmaxd, irmind:irm), intent (out) :: vnspll
    ! .. Local Scalars
    integer :: ir, j, lm1, lm2, lm3
    ! ..
    vnspll(1:lmmaxd, 1:lmmaxd, irmind:irm) = 0.0e0_dp

    do j = 1, iend
      lm1 = icleb(j, 1)
      lm2 = icleb(j, 2)
      lm3 = icleb(j, 3)
      !> @warning In the old version of KKRimp the following is only done for lm1<=lmmaxd && lm2<=lmmaxd && lm3>1 @endwarning
      do ir = irmind, irm
        vnspll(lm1, lm2, ir) = vnspll(lm1, lm2, ir) + cleb(j, 1)*vins(ir, lm3)
      end do
    end do
    ! ----------------------------------------------------------------------------
    ! Use symmetry of the gaunt coef.
    ! ----------------------------------------------------------------------------
    do lm1 = 1, lmmaxd
      do lm2 = 1, lm1 - 1
        do ir = irmind, irm
          vnspll(lm2, lm1, ir) = vnspll(lm1, lm2, ir)
        end do
      end do
    end do

    !> @warning The following part is commented out in the old version of KKRimp routine! @endwarning
    do lm1 = 1, lmmaxd
      do ir = irmind, irm
        vnspll(lm1, lm1, ir) = vnspll(lm1, lm1, ir) + vins(ir, 1)
      end do
    end do

  end subroutine vllns

end module mod_vllns
