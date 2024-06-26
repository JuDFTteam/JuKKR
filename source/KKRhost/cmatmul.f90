!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_cmatmul
  
  private
  public :: cmatmul

contains

  !-------------------------------------------------------------------------------
  !> Summary: Complex matrix-matrix multiplication
  !> Author: 
  !> Category: KKRhost, undefined
  !> Deprecated: True ! This needs to be set to True for deprecated subroutines
  !>
  !> Perform  the matrix-matrix operation           C = A * B
  !>
  !>  A,B,C   complex  SQUARE  N x N - matrices
  !>  N       dimension of A, B and C
  !>  M       array size of A, B, C with M >= N
  !>
  !> @note
  !> Can easily be replaced with intrinsic MATMUL or LAPACK call (should be much more optimized)
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine cmatmul(n, m, a, b, c)
    use :: mod_datatypes, only: dp
    use :: mod_constants, only: czero
    implicit none

    ! PARAMETER definitions
    real (kind=dp), parameter :: eps = 1.0e-12_dp

    ! Dummy arguments
    integer :: m, n
    complex (kind=dp) :: a(m, m), b(m, m), c(m, m)

    ! Local variables
    complex (kind=dp) :: blj
    integer :: i, j, l

    do j = 1, n
      do i = 1, n
        c(i, j) = czero
      end do
    end do

    do j = 1, n
      do l = 1, n
        blj = b(l, j)
        if (abs(blj)>eps) then
          do i = 1, n
            c(i, j) = c(i, j) + a(i, l)*blj
          end do
        end if
      end do
    end do

  end subroutine cmatmul

end module mod_cmatmul
