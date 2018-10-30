      SUBROUTINE RINTSIMP(FX,JTOP,CINT)
C   ********************************************************************
C   *                                                                  *
C   *  SIMPSON - INTERGRATION FOR  REAL   INTEGRAND  FX FROM 1 TO JTOP *
C   *  AND EQUIDISTANT MESH    I                                       *
C   *   INT = [ F1 + 4*F2 + 2*F3 + .... + 4F(N-1) + FN ]/3     N:ODD   *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C
C Dummy arguments
C
      REAL*8 CINT
      INTEGER JTOP
      REAL*8 FX(JTOP)
C
C Local variables
C
      INTEGER I
      REAL*8 SIMP
C
      CINT = FX(1)
      SIMP = -1.0D0
C
      DO I = 2,JTOP - 1
         SIMP = -SIMP
         CINT = CINT + (3.0D0+SIMP)*FX(I)
      END DO
C
      CINT = (CINT+FX(JTOP))/3.0D0
C
      END
