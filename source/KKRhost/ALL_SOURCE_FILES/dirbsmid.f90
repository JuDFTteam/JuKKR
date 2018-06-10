subroutine dirbsmid(y, dydx, nv, xs, htot, nstep, yout, b, v, r, drdi, nmesh)
!   ********************************************************************
!   *                                                                  *
!   *   modified midpoint step to support the  Burlisch-Stoer method   *
!   *   on exit:  the incremented variable is in   YOUT                *
!   *                                                                  *
!   *   see: numerical recipes chapter 15.3                            *
!   *                                                                  *
!   ********************************************************************

  implicit none

  include 'sprkkr_rmesh.dim'

! Dummy arguments
  real *8 :: htot, xs
  integer :: nmesh, nstep, nv
  real *8 :: b(nrmax), drdi(nrmax), r(nrmax), v(nrmax)
  complex *16 :: dydx(nv), y(nv), yout(nv)

! Local variables
  real *8 :: h, h2, x
  integer :: i, n
  complex *16 :: swap, ym(ncfmax), yn(ncfmax)

  h = htot/nstep
  do i = 1, nv
    ym(i) = y(i)
    yn(i) = y(i) + h*dydx(i)
  end do
  x = xs + h

  call dirbsrad(x, yn, yout, drdi, b, v, r, nmesh)

  h2 = 2.d0*h
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
    yout(i) = 0.5d0*(ym(i)+yn(i)+h*yout(i))

  end do

end subroutine
