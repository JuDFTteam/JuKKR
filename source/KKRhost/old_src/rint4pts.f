      SUBROUTINE RINT4PTS(Y,JTOP,Z)
C
C   ********************************************************************
C   *                                                                  *
C   *      perform the integral  Z(i)   =  INT   Y(i') di'             *
C   *                                    R=0..R(i)                     *
C   *                                                                  *
C   *      via a 4-point integration formula                           *
C   *                                                                  *
C   *      JTOP:     Y is tabulated form 1 .. JTOP                     *
C   *      Y(i):     function to be integrated                         *
C   *                                                                  *
C   *                       REAL    - VERSION                          *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE

C
C Dummy arguments
C
      INTEGER JTOP
      REAL*8 Y(JTOP),Z(JTOP)
C
C Local variables
C
      INTEGER I,IG,J,K,M,N1,N2
      REAL*8 Q(5,5),Q5(5,5),S,SVN
C
      DATA Q5/0.D0,251.D0,232.D0,243.D0,224.D0,0.D0,646.D0,992.D0,
     &     918.D0,1024.D0,0.D0, - 264.D0,192.D0,648.D0,384.D0,0.D0,
     &     106.D0,32.D0,378.D0,1024.D0,0.D0, - 19.D0, - 8.D0, - 27.D0,
     &     224.D0/
C
      DO I = 1,5
         DO J = 1,5
            Q(I,J) = Q5(I,J)/720.0D0
         END DO
      END DO
C
      Z(1) = 0.0D0
      SVN = Z(1)
C
      DO IG = 1,JTOP - 4,4
         N1 = IG
         N2 = IG + 4
         DO M = N1 + 1,N2
            I = M - N1 + 1
            S = SVN
            DO K = N1,N2
               J = K - N1 + 1
               S = S + Q(I,J)*Y(K)
            END DO
            Z(M) = S
         END DO
         SVN = Z(N2)
      END DO
C
      IF ( N2.NE.JTOP ) THEN
         N1 = JTOP - 4
         N2 = JTOP
         SVN = Z(N1)
         DO M = N1 + 1,N2
            I = M - N1 + 1
            S = SVN
            DO K = N1,N2
               J = K - N1 + 1
               S = S + Q(I,J)*Y(K)
            END DO
            Z(M) = S
         END DO
      END IF
C
      END
