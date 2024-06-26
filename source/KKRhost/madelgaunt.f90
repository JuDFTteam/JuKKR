!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Calculation of the Gaunt coefficients
!> Author: 
!> Calculation of the Gaunt coefficients
!------------------------------------------------------------------------------------
module mod_madelgaunt
  use :: mod_datatypes, only: dp
  private :: dp

contains
  !-------------------------------------------------------------------------------
  !> Summary: Calculation of the Gaunt coefficients
  !> Author: 
  !> Category: electrostatics, KKRhost 
  !> Deprecated: False
  !> Calculation of the Gaunt coefficients
  !-------------------------------------------------------------------------------
  !> @note Attention: Dimension NCLEBD appears sometimes as `NCLEB1` an empirical 
  !> factor - it has to be optimized. 
  !> 
  !> Jonathan Chico 21.09.2018: Unsure if previous note is still valid
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine madelgaunt(lpot,yrg,wg,cleb,icleb,iend,lassld,nclebd)
    implicit none
    real (kind=dp), parameter :: eps = 1.0e-12_dp
    ! ..
    ! .. Input variables
    integer, intent(in) :: lpot   !! Maximum l component in potential expansion
    integer, intent(in) :: lassld !! 4*lmax
    integer, intent(in) :: nclebd !! (LMAX*2+1)**2 * (LMAX+1)**2
    real (kind=dp), dimension(lassld), intent(in) :: wg !! Integr. weights for Legendre polynomials
    real (kind=dp), dimension(lassld, 0:lassld, 0:lassld), intent(in) :: yrg !! Spherical harmonics (GAUNT2)
    ! ..
    ! .. Output variables 
    integer, intent(out) :: iend  !! Number of nonzero gaunt coefficients
    integer, dimension(nclebd,3), intent(out) :: icleb !! Pointer array
    real (kind=dp), dimension(nclebd), intent(out) :: cleb !! GAUNT coefficients (GAUNT)
    ! ..
    ! .. Local scalars
    real (kind=dp) :: clecg, factor, s
    integer :: i, j, l1, l2, l3, m1, m1a, m1s, m2, m2a, m2s, m3, m3a, m3s
    ! ..
    ! .. Intrinsic functions
    intrinsic :: abs, real, sign

    ! --> set up of the gaunt coefficients with an index field
    ! recognize that they are needed here only for l3=l1+l2

    if (2*lpot>lassld) then
      write (6, *) 'Dim ERROR in MADELGAUNT -- 2*LPOT > LASSLD', 2*lpot, lassld
      stop
    end if

    i = 1
    ! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
    do l1 = 0, lpot
      do l2 = 0, lpot
        l3 = l1 + l2
        ! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
        do m1 = -l1, l1
          do m2 = -l2, l2
            do m3 = -l3, l3
              m1s = sign(1, m1)
              m2s = sign(1, m2)
              m3s = sign(1, m3)
              ! **********************************************************************
              if (m1s*m2s*m3s>=0) then
                m1a = abs(m1)
                m2a = abs(m2)
                m3a = abs(m3)

                factor = 0.0e0_dp
                if (m1a+m2a==m3a) factor = factor + real(3*m3s+sign(1,-m3), kind=dp)/8.0e0_dp
                if (m1a-m2a==m3a) factor = factor + real(m1s, kind=dp)/4.0e0_dp
                if (m2a-m1a==m3a) factor = factor + real(m2s, kind=dp)/4.0e0_dp
                ! ======================================================================
                if (abs(factor)>eps) then
                  if (m1s*m2s/=1 .or. m2s*m3s/=1 .or. m1s*m3s/=1) factor = -factor

                  s = 0.0e0_dp
                  do j = 1, lassld
                    s = s + wg(j)*yrg(j, l1, m1a)*yrg(j, l2, m2a)*yrg(j, l3, m3a)
                  end do

                  clecg = s*factor
                  ! ----------------------------------------------------------------------
                  if (abs(clecg)>1.e-10_dp) then
                    cleb(i) = clecg
                    icleb(i, 1) = l1*(l1+1) + m1 + 1
                    icleb(i, 2) = l2*(l2+1) + m2 + 1
                    icleb(i, 3) = l3*(l3+1) + m3 + 1
                    i = i + 1
                    if (i>nclebd) then
                      write (6, fmt='(2I10)') i, nclebd
                      stop ' Dim stop in MADELGAUNT '
                    end if
                  end if
                  ! ----------------------------------------------------------------------
                end if
                ! ======================================================================
              end if
              ! **********************************************************************
            end do
          end do
        end do
        ! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
      end do
    end do
    ! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
    iend = i - 1
  end subroutine madelgaunt

end module mod_madelgaunt