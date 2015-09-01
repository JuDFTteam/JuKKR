C*==checkrmat.f    processed by SPAG 6.05Rc at 12:32 on  3 Oct 2021
      LOGICAL FUNCTION CHECKRMAT(RMAT,CO1,SI1,CO2,SI2,CO3,SI3,I,J)
C   ********************************************************************
C   *                                                                  *
C   *  check whether the values of the cosinus and sinus found for the *
C   *  Euler angles TET1, TET2, TET3 are consistent with the           *
C   *  rotation matrix   RMAT                                          *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      DOUBLE PRECISION CO1,CO2,CO3,SI1,SI2,SI3
      INTEGER I,J
      DOUBLE PRECISION RMAT(3,3)
C
C Local variables
C
      REAL*8 A,B
      LOGICAL EQUAL
C
C*** End of declarations rewritten by SPAG
C
      EQUAL(A,B) = (ABS(A-B).LT.1D-7)
C
      CHECKRMAT = .FALSE.
C
      IF ( I.EQ.1 ) THEN
         IF ( J.EQ.1 ) THEN
            CHECKRMAT = EQUAL(RMAT(1,1),CO3*CO2*CO1-SI3*SI1)
         ELSE IF ( J.EQ.2 ) THEN
            CHECKRMAT = EQUAL(RMAT(1,2),CO3*CO2*SI1+SI3*CO1)
         ELSE IF ( J.EQ.3 ) THEN
            CHECKRMAT = EQUAL(RMAT(1,3),-CO3*SI2)
         END IF
      ELSE IF ( I.EQ.2 ) THEN
         IF ( J.EQ.1 ) THEN
            CHECKRMAT = EQUAL(RMAT(2,1),-SI3*CO2*CO1-CO3*SI1)
         ELSE IF ( J.EQ.2 ) THEN
            CHECKRMAT = EQUAL(RMAT(2,2),-SI3*CO2*SI1+CO3*CO1)
         ELSE IF ( J.EQ.3 ) THEN
            CHECKRMAT = EQUAL(RMAT(2,3),SI3*SI2)
         END IF
      ELSE IF ( J.EQ.1 ) THEN
         CHECKRMAT = EQUAL(RMAT(3,1),SI2*CO1)
      ELSE IF ( J.EQ.2 ) THEN
         CHECKRMAT = EQUAL(RMAT(3,2),SI2*SI1)
      ELSE IF ( J.EQ.3 ) THEN
         CHECKRMAT = EQUAL(RMAT(3,3),CO2)
      END IF
      END
