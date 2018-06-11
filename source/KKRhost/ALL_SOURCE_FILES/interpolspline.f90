    Subroutine interpolspline(rmesh, rmeshnew, vpot, vpotnew, nrmax, nrmaxnew)
      Use mod_datatypes, Only: dp
      Implicit None
!interface
      Integer :: nrmax
      Integer :: nrmaxnew
      Real (Kind=dp) :: rmesh(nrmax)
      Real (Kind=dp) :: rmeshnew(nrmaxnew)
      Real (Kind=dp) :: vpot(nrmax)
      Real (Kind=dp) :: vpotnew(nrmaxnew)
!local
      Real (Kind=dp) :: maxa
      Real (Kind=dp) :: spline(nrmax)
      Real (Kind=dp) :: parsum, parsumderiv, r0
      Integer :: ir

      maxa = 1.E35_dp
      Call spline_real(nrmax, rmesh, vpot, nrmax, maxa, maxa, spline)
!           CALL SPLINE(IRMDJJ,R,VM2Z,NR,maxa,maxa,VM2ZB)

      Do ir = 1, nrmaxnew
        r0 = rmeshnew(ir)
        Call splint_real(rmesh, vpot, spline, nrmax, r0, parsum, parsumderiv)
        vpotnew(ir) = parsum
      End Do
    End Subroutine
