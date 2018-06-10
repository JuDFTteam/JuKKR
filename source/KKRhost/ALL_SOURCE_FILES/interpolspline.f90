subroutine interpolspline(rmesh, rmeshnew, vpot, vpotnew, nrmax, nrmaxnew)
  implicit none
!interface
  integer :: nrmax
  integer :: nrmaxnew
  double precision :: rmesh(nrmax)
  double precision :: rmeshnew(nrmaxnew)
  double precision :: vpot(nrmax)
  double precision :: vpotnew(nrmaxnew)
!local
  double precision :: maxa
  double precision :: spline(nrmax)
  double precision :: parsum, parsumderiv, r0
  integer :: ir

  maxa = 1.d35
  call spline_real(nrmax, rmesh, vpot, nrmax, maxa, maxa, spline)
!           CALL SPLINE(IRMDJJ,R,VM2Z,NR,maxa,maxa,VM2ZB)

  do ir = 1, nrmaxnew
    r0 = rmeshnew(ir)
    call splint_real(rmesh, vpot, spline, nrmax, r0, parsum, parsumderiv)
    vpotnew(ir) = parsum
  end do
end subroutine

