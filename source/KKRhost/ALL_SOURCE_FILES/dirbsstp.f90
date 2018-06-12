subroutine dirbsstp(y, dydx, nv, x, htry, eps, yscal, b, v, r, drdi, nmesh)
  ! ********************************************************************
  ! *                                                                  *
  ! *   Burlisch-Stoer step with monitoring of local truncation error  *
  ! *   on entry: X,Y,DXDY  for last mesh-point                        *
  ! *   on exit:  X,Y,DXDY  updated for X = X(last) + HTRY             *
  ! *                                                                  *
  ! *   see: numerical recipes chapter 15.4                            *
  ! *                                                                  *
  ! *   note: don't set NUSE    > NUSEMAX in <DIRBSRZE>                *
  ! *         don't set ISEQMAX > ISEQMAX in <DIRBSRZE>                *
  ! *         no step size adjusted in case of no convergency > STOP   *
  ! *                                                                  *
  ! ********************************************************************

  use :: mod_datatypes, only: dp
  use :: mod_types, only: t_inc
  implicit none

  include 'sprkkr_rmesh.dim'

  ! PARAMETER definitions
  integer :: iseqmax, nuse
  parameter (iseqmax=30, nuse=7)
  complex (kind=dp) :: tiny
  ! Bereshad:
  ! parameter (  tiny  = (1.0d-20,1.0d-20)  )
  ! dadurch wird der Imaginaerteil fuer die Skalierungsfunktion nicht 00
  ! klein.
  ! Ich hatte da spruenge in dem errmax; das verlaeuft jetzt stetig, hat aber
  ! mit
  ! den kleinen Werten doch so seine Konvergenzprobleme.

  parameter (tiny=(1.0d-20,1.0d-20))

  ! Dummy arguments
  real (kind=dp) :: eps, htry, x
  integer :: nmesh, nv
  real (kind=dp) :: b(nrmax), drdi(nrmax), r(nrmax), v(nrmax)
  complex (kind=dp) :: dydx(nv), y(nv), yscal(nv)

  ! Local variables
  complex (kind=dp) :: dysav(ncfmax), yerr(ncfmax), ysav(ncfmax), yseq(ncfmax)
  real (kind=dp) :: errmax, h, xest, xsav
  integer :: i, j, nseq(iseqmax)

  data nseq/2, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, &
    768, 1024, 1536, 2048, 3072, 4096, 6144, 8192, 12288, 16384, 24576, 32768, &
    49152, 65536/

  h = htry
  xsav = x

  do i = 1, nv
    ysav(i) = y(i)
    dysav(i) = dydx(i)
  end do
  do i = 1, iseqmax

    call dirbsmid(ysav, dysav, nv, xsav, h, nseq(i), yseq, b, v, r, drdi, &
      nmesh)
    xest = (h/nseq(i))**2

    call dirbsrze(i, xest, yseq, y, yerr, nv, nuse)

    errmax = 0.0d0
    do j = 1, nv
      errmax = dble(max(errmax,abs(yerr(j)/(yscal(j)+tiny))))
    end do
    errmax = errmax/eps
    if (errmax<1.0d0) then
      x = x + h

      call dirbsrad(x, y, dydx, drdi, b, v, r, nmesh)


      return
    end if
  end do

  if (errmax<1000d0) then
    x = x + h
    call dirbsrad(x, y, dydx, drdi, b, v, r, nmesh)
    if (t_inc%i_write>0) then
      write (1337, *) '<DIRBSSTP>  not converged after ', iseqmax, &
        ' refinements'
      write (1337, *) 'step size will not be adjusted !!!!!!'
      write (1337, *) 'max. relative error : ', errmax*eps
      write (1337, *) 'tolerance             ', eps
      write (1337, *) 'grid position  X      ', x
    end if
  else
    if (t_inc%i_write>0) then
      write (6, *) '<DIRBSSTP>  not converged after ', iseqmax, ' refinements'
      write (6, *) 'step size will not be adjusted !!!!!!'
      write (6, *) 'max. relative error : ', errmax*eps
      write (6, *) 'tolerance             ', eps
      write (6, *) 'grid position  X      ', x
    end if
    stop
  end if


  return
end subroutine dirbsstp
