!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Invert a matrix \(A\) using the Gauss-Jordan algorithm.
!> Author: 
!> Invert a matrix \(A\) using the Gauss-Jordan algorithm. The 1- matrix is not 
!> set up and use is made of its structure.
!------------------------------------------------------------------------------------
module mod_rinvgj
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Invert a matrix \(A\) using the Gauss-Jordan algorithm.
  !> Author: 
  !> Category: numerical-tools, KKRhost
  !> Deprecated: False 
  !> Invert a matrix \(A\) using the Gauss-Jordan algorithm. The 1- matrix is not 
  !> set up and use is made of its structure.
  !-------------------------------------------------------------------------------
  subroutine rinvgj(ainv, a, arraydim, n)

    implicit none

    ! Dummy arguments
    integer, intent(in) :: n !! Number of columns
    integer, intent(in) :: arraydim !! Dimension of the array
    real (kind=dp), dimension(arraydim, arraydim), intent(inout) :: a   !! Matrix to be inverted
    real (kind=dp), dimension(arraydim, arraydim), intent(out) :: ainv  !! Inverted matrix

    ! Local variables
    integer :: icol, l, ll
    real (kind=dp) :: t, t1

    ainv(1, 1) = 0e0_dp
    ! scan columns
    do icol = 1, n

      ! make A(ICOL,ICOL) = 1
      t1 = 1.0e0_dp/a(icol, icol)
      do l = (icol+1), n
        a(icol, l) = a(icol, l)*t1
      end do

      do l = 1, (icol-1)
        ainv(icol, l) = ainv(icol, l)*t1
      end do
      ainv(icol, icol) = t1

      ! make A(LL,ICOL) = 0 for LL<>ICOL
      do ll = 1, n
        if (ll/=icol) then
          t = a(ll, icol)
          do l = (icol+1), n
            a(ll, l) = a(ll, l) - a(icol, l)*t
          end do

          do l = 1, (icol-1)
            ainv(ll, l) = ainv(ll, l) - ainv(icol, l)*t
          end do
          ainv(ll, icol) = -t1*t
        end if
      end do
    end do

  end subroutine rinvgj

end module mod_rinvgj
