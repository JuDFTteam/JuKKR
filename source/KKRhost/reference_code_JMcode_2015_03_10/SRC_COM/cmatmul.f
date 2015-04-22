      SUBROUTINE CMATMUL(N,M,A,B,C)
C   ********************************************************************
C   *                                                                  *
C   *   perform  the matrix-matrix operation           C = A * B       *
C   *                                                                  *
C   *   A,B,C   complex  SQUARE  N x N - matrices                      *
C   *   N       dimension of A, B and C                                *
C   *   M       array size of A, B, C with M >= N                      *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT DOUBLE COMPLEX(A-H,O-Z)
C
C PARAMETER definitions
C
      DOUBLE COMPLEX C0
      PARAMETER (C0=(0.0D0,0.0D0))
C
C Dummy arguments
C
      INTEGER M,N
      DOUBLE COMPLEX A(M,M),B(M,M),C(M,M)
C
C Local variables
C
      DOUBLE COMPLEX BLJ
      INTEGER I,J,L
C
      DO J = 1,N
         DO I = 1,N
            C(I,J) = C0
         END DO
      END DO
C
      DO J = 1,N
         DO L = 1,N
            BLJ = B(L,J)
            IF ( BLJ.NE.C0 ) THEN
               DO I = 1,N
                  C(I,J) = C(I,J) + A(I,L)*BLJ
               END DO
            END IF
         END DO
      END DO
C
      END
