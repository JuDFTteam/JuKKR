SUBROUTINE interpolspline(rmesh,rmeshnew,vpot,vpotnew,  &
        nrmax,nrmaxnew)
IMPLICIT NONE
!interface
INTEGER :: nrmax
INTEGER :: nrmaxnew
DOUBLE PRECISION :: rmesh(nrmax)
DOUBLE PRECISION :: rmeshnew(nrmaxnew)
DOUBLE PRECISION :: vpot(nrmax)
DOUBLE PRECISION :: vpotnew(nrmaxnew)
!local
DOUBLE PRECISION :: maxa
DOUBLE PRECISION :: spline(nrmax)
DOUBLE PRECISION :: parsum, parsumderiv,r0
INTEGER :: ir
maxa = 1.d35
CALL spline_real(nrmax,rmesh,vpot,nrmax,maxa,maxa,spline)
!           CALL SPLINE(IRMDJJ,R,VM2Z,NR,maxa,maxa,VM2ZB)

DO ir = 1,nrmaxnew
  r0 = rmeshnew(ir)
  CALL splint_real(rmesh,vpot,spline,nrmax,r0,parsum,parsumderiv)
  vpotnew(ir) = parsum
END DO
END SUBROUTINE  interpolspline

