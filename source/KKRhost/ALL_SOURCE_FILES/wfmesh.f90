SUBROUTINE wfmesh(e,ek,cvlight,nsra,z,r,s,rs,irm,irmd,lmaxd)
!.. Scalar Arguments ..
      DOUBLE COMPLEX E,EK
      DOUBLE PRECISION CVLIGHT,Z
      INTEGER IRM,IRMD,LMAXD,NSRA
!..
!.. Intrinsic Functions ..
      INTRINSIC DBLE,SQRT
!..
!.. Array Arguments ..
      DOUBLE PRECISION R(IRMD),RS(IRMD,0:LMAXD),S(0:LMAXD)
!..
!.. Local Scalars ..
      DOUBLE PRECISION S1
      INTEGER IR,L
!..
IF (nsra == 1) ek = SQRT(e)
IF (nsra == 2) ek = SQRT(e+e*e/ (cvlight*cvlight))
DO l = 0,lmaxd
  
  IF (nsra == 2) THEN
    s1 = SQRT(DBLE(l*l+l+1)-4.0D0*z*z/ (cvlight*cvlight))
    IF (z == 0.0D0) s1 = DBLE(l)
  ELSE
    s1 = DBLE(l)
  END IF
  s(l) = s1
  rs(1,l) = 0.0D0
  DO ir = 2,irm
    rs(ir,l) = r(ir)**s1
  END DO
  DO ir = irm+1, irmd
    rs(ir,l) = 0.0D0
  END DO
  
END DO
RETURN
END SUBROUTINE wfmesh
