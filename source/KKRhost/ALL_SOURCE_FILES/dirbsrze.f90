SUBROUTINE dirbsrze(iest,xest,yest,yz,dy,nv,nuse)
!   ********************************************************************
!   *                                                                  *
!   *   diagonal rational function extrapolation to support the        *
!   *   Burlisch-Stoer method                                          *
!   *                                                                  *
!   *   see: numerical recipes chapter 15.4                            *
!   *                                                                  *
!   ********************************************************************

IMPLICIT NONE

! PARAMETER definitions
INTEGER NCFMAX,ISEQMAX,NUSEMAX
PARAMETER (NCFMAX=8,ISEQMAX=30,NUSEMAX=7)

! Dummy arguments
INTEGER IEST,NUSE,NV
REAL*8 XEST
COMPLEX*16 DY(NV),YEST(NV),YZ(NV)

! Local variables
COMPLEX*16 B,B1,C,D(NCFMAX,NUSEMAX),DDY,V,YY
REAL*8 FX(NUSEMAX),X(ISEQMAX)
INTEGER J,K,M1
SAVE B,B1,C,D,DDY,FX,J,K,M1,V,X,YY

x(iest) = xest
IF ( iest == 1 ) THEN
  DO j = 1,nv
    yz(j) = yest(j)
    d(j,1) = yest(j)
    dy(j) = yest(j)
  END DO
ELSE
  m1 = MIN(iest,nuse)
  DO k = 1,m1 - 1
    fx(k+1) = x(iest-k)/xest
  END DO
  DO j = 1,nv
    yy = yest(j)
    v = d(j,1)
    c = yy
    d(j,1) = yy
    DO k = 2,m1
      b1 = fx(k)*v
      b = b1 - c
      IF ( b /= 0. ) THEN
        b = (c-v)/b
        ddy = c*b
        c = b1*b
      ELSE
        ddy = v
      END IF
      v = d(j,k)
      d(j,k) = ddy
      yy = yy + ddy
    END DO
    dy(j) = ddy
    yz(j) = yy
  END DO
END IF
END SUBROUTINE dirbsrze
