!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_dscal1
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Scaling of vector by constant
  !> Author: 
  !> Category: KKRhost, numerical-tools
  !> Deprecated: True ! This needs to be set to True for deprecated subroutines
  !>
  !> @note
  !> This can be easily be replace with `dx(:) = a*dx(:)` but this routine is anyways only used by
  !> dead parts of the core base
  !> @endnote
  !>
  !>  Scales a vector by a constant  dx(i) -> a * dx(i)
  !> ----------------------------------------------------------------------
  !> i Inputs:
  !> i   n     :lenght of dx and dy
  !> i   da    :constant
  !> i   dx    :vector
  !> i   incx  :incrementation for x
  !> o Outputs:
  !> o   dx    :vector
  !> r Remarks:
  !> r   Adapted from: jack dongarra, linpack, 3/11/78.
  !-------------------------------------------------------------------------------
  subroutine dscal1(n, da, dx, incx)

    implicit none
    ! Passed parameters:
    real (kind=dp) :: da, dx(*)
    integer :: incx, n
    ! Local parameters:
    integer :: i, m, mp1, nincx

    if (n<=0 .or. incx<=0) return
    if (incx/=1) then
      ! ----- code for increment not equal to 1
      nincx = n*incx
      do i = 1, nincx, incx
        dx(i) = da*dx(i)
      end do
    else
      ! ----- code for increment equal to 1
      m = mod(n, 5)
      if (m/=0) then
        do i = 1, m
          dx(i) = da*dx(i)
        end do
        if (n<5) return
      end if
      mp1 = m + 1
      do i = mp1, n, 5
        dx(i) = da*dx(i)
        dx(i+1) = da*dx(i+1)
        dx(i+2) = da*dx(i+2)
        dx(i+3) = da*dx(i+3)
        dx(i+4) = da*dx(i+4)
      end do
    end if
  end subroutine dscal1

end module mod_dscal1
