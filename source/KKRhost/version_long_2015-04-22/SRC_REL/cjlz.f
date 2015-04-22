      FUNCTION CJLZ(L,Z)
C   ********************************************************************
C   *                                                                  *
C   *   SPHERICAL BESSEL-FUNCTION  J(L,Z)  FOR COMPLEX ARGUMENT  Z     *
C   *                  see:  e.g. MERZBACHER EQ. (10.22)               *
C   *                                                                  *
C   ********************************************************************
C

      IMPLICIT NONE

C
C PARAMETER definitions
C
      COMPLEX*16 C1
      PARAMETER (C1=(1.0D0,0.0D0))
      INTEGER LP2MAX
      PARAMETER (LP2MAX=25)

C
C Dummy arguments
C
      INTEGER L
      COMPLEX*16 Z
      COMPLEX*16 CJLZ
C
C Local variables
C
      
      REAL*8 DFAC
      COMPLEX*16 DT,S(LP2MAX),T,ZSQ
      INTEGER I,K,LLP1
C
      ZSQ = Z*Z
      LLP1 = L + L + 1
C
      IF ( ABS(ZSQ/DBLE(LLP1)).LE.10.D0 ) THEN
C
         DFAC = 1.0D0
         DO K = 3,LLP1,2
            DFAC = DFAC*DBLE(K)
         END DO
C
         DT = C1
         T = C1
         DO I = 2,400,2
            DT = -DT*ZSQ/DBLE(I*(I+LLP1))
            T = T + DT
            IF ( ABS(DT).LT.1.0D-10 ) GOTO 50
         END DO
C
 50      CONTINUE
         CJLZ = T*Z**L/DFAC
C
      ELSE
         IF ( L.GT.23 ) STOP '<cjlz>: l too large'
C
         S(2) = SIN(Z)/Z
         IF ( L.LE.0 ) THEN
            CJLZ = S(2)*Z**L
            RETURN
         END IF
C
         S(1) = COS(Z)
         DO I = 3,L + 2
            S(I) = (S(I-1)*(2*I-5)-S(I-2))/ZSQ
         END DO
         CJLZ = S(L+2)*Z**L
C
      END IF
      END 
