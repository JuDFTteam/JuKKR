    Subroutine dirbsmid(y, dydx, nv, xs, htot, nstep, yout, b, v, r, drdi, &
      nmesh)
      Use mod_datatypes, Only: dp
!   ********************************************************************
!   *                                                                  *
!   *   modified midpoint step to support the  Burlisch-Stoer method   *
!   *   on exit:  the incremented variable is in   YOUT                *
!   *                                                                  *
!   *   see: numerical recipes chapter 15.3                            *
!   *                                                                  *
!   ********************************************************************

      Implicit None

      Include 'sprkkr_rmesh.dim'

! Dummy arguments
      Real (Kind=dp) :: htot, xs
      Integer :: nmesh, nstep, nv
      Real (Kind=dp) :: b(nrmax), drdi(nrmax), r(nrmax), v(nrmax)
      Complex (Kind=dp) :: dydx(nv), y(nv), yout(nv)

! Local variables
      Real (Kind=dp) :: h, h2, x
      Integer :: i, n
      Complex (Kind=dp) :: swap, ym(ncfmax), yn(ncfmax)

      h = htot/nstep
      Do i = 1, nv
        ym(i) = y(i)
        yn(i) = y(i) + h*dydx(i)
      End Do
      x = xs + h

      Call dirbsrad(x, yn, yout, drdi, b, v, r, nmesh)

      h2 = 2.E0_dp*h
      Do n = 2, nstep
        Do i = 1, nv
          swap = ym(i) + h2*yout(i)
          ym(i) = yn(i)
          yn(i) = swap
        End Do
        x = x + h
        Call dirbsrad(x, yn, yout, drdi, b, v, r, nmesh)
      End Do
      Do i = 1, nv
        yout(i) = 0.5E0_dp*(ym(i)+yn(i)+h*yout(i))

      End Do

    End Subroutine
