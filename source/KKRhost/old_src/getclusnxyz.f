C*==getclusnxyz.f    processed by SPAG 6.05Rc at 15:44 on 18 Oct 2004
      SUBROUTINE GETCLUSNXYZ(CLURAD,BRAVAIS,NDIM,CLURADSQ,NBR)
C **********************************************************************
C *                                                                    *
C * Given a spherical cluster of radius CLURAD it determines the three *
C * integers N1,N2,N3 such that any vector                             *
C *                                                                    *
C *    R_i = r_i + SUM_j  N_j * a_j                                    *
C *                                                                    *
C *  with i = 1,NAEZ and a_j the primitive Bravais vectors, is inside  *
C *  the cluster. Subroutine also returns the CLURAD**2 value          *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
C ..
C ..  Arguments
      DOUBLE PRECISION CLURAD,CLURADSQ
      DOUBLE PRECISION BRAVAIS(3,3)
      INTEGER NDIM,NBR(3)
C .. 
C ..  Locals
      DOUBLE PRECISION DR(3)
      INTEGER I,J
      INTEGER INT
C ..
      DO I = 1,NDIM
         DR(I) = 0D0
         DO J = 1,NDIM
            DR(I) = DR(I) + BRAVAIS(J,I)*BRAVAIS(J,I)
         END DO
         DR(I) = SQRT(DR(I))
      END DO
C
      IF ( ABS(CLURAD).LT.1D-6 ) THEN
         DO I = 1,NDIM
            NBR(I) = 0
         END DO
         CLURADSQ = 1D10
      ELSE
         DO I = 1,NDIM
            NBR(I) = INT(CLURAD/DR(I)) + 2
         END DO
         CLURADSQ = CLURAD*CLURAD
      END IF
C
      END
