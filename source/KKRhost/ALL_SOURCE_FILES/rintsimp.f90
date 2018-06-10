SUBROUTINE rintsimp(fx,jtop,cint)
!   ********************************************************************
!   *                                                                  *
!   *  SIMPSON - INTERGRATION FOR  REAL   INTEGRAND  FX FROM 1 TO JTOP *
!   *  AND EQUIDISTANT MESH    I                                       *
!   *   INT = [ F1 + 4*F2 + 2*F3 + .... + 4F(N-1) + FN ]/3     N:ODD   *
!   *                                                                  *
!   ********************************************************************

IMPLICIT NONE

! Dummy arguments
REAL*8 CINT
INTEGER JTOP
REAL*8 FX(JTOP)

! Local variables
INTEGER I
REAL*8 SIMP

cint = fx(1)
simp = -1.0D0

DO i = 2,jtop - 1
  simp = -simp
  cint = cint + (3.0D0+simp)*fx(i)
END DO

cint = (cint+fx(jtop))/3.0D0

END SUBROUTINE rintsimp
