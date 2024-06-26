!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_dirbsmid

contains

  !-------------------------------------------------------------------------------
  !> Summary: Burlisch-Stoer method
  !> Author: 
  !> Category: KKRhost, dirac, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Modified midpoint step to support the Burlisch-Stoer method
  !> on exit:  the incremented variable is in   YOUT
  !>
  !> see: numerical recipes chapter 15.3
  !-------------------------------------------------------------------------------
  subroutine dirbsmid(y, dydx, nv, xs, htot, nstep, yout, b, v, r, drdi, nmesh)

    use :: mod_datatypes, only: dp
    use :: mod_dirbsrad, only: dirbsrad
    implicit none

    include 'sprkkr_rmesh.dim'

    ! Dummy arguments
    real (kind=dp) :: htot, xs
    integer :: nmesh, nstep, nv
    real (kind=dp) :: b(nrmax), drdi(nrmax), r(nrmax), v(nrmax)
    complex (kind=dp) :: dydx(nv), y(nv), yout(nv)

    ! Local variables
    real (kind=dp) :: h, h2, x
    integer :: i, n
    complex (kind=dp) :: swap, ym(ncfmax), yn(ncfmax)

    h = htot/nstep
    do i = 1, nv
      ym(i) = y(i)
      yn(i) = y(i) + h*dydx(i)
    end do
    x = xs + h

    call dirbsrad(x, yn, yout, drdi, b, v, r, nmesh)

    h2 = 2.e0_dp*h
    do n = 2, nstep
      do i = 1, nv
        swap = ym(i) + h2*yout(i)
        ym(i) = yn(i)
        yn(i) = swap
      end do
      x = x + h
      call dirbsrad(x, yn, yout, drdi, b, v, r, nmesh)
    end do
    do i = 1, nv
      yout(i) = 0.5e0_dp*(ym(i)+yn(i)+h*yout(i))

    end do

  end subroutine dirbsmid

end module mod_dirbsmid
