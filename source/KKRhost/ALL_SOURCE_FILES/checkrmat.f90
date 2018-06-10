FUNCTION checkrmat(rmat,co1,si1,co2,si2,co3,si3,i,j)
!   ********************************************************************
!   *                                                                  *
!   *  check whether the values of the cosinus and sinus found for the *
!   *  Euler angles TET1, TET2, TET3 are consistent with the           *
!   *  rotation matrix   RMAT                                          *
!   *                                                                  *
!   ********************************************************************

IMPLICIT NONE

! Dummy arguments
DOUBLE PRECISION CO1,CO2,CO3,SI1,SI2,SI3
INTEGER I,J
LOGICAL CHECKRMAT
DOUBLE PRECISION RMAT(3,3)

! Local variables
DOUBLE PRECISION A,B
LOGICAL EQUAL
LOGICAL RESULT

equal(a,b) = (ABS(a-b) < 1D-7)

result = .false.

IF ( i == 1 ) THEN
  IF ( j == 1 ) THEN
    result = equal(rmat(1,1),co3*co2*co1-si3*si1)
  ELSE IF ( j == 2 ) THEN
    result = equal(rmat(1,2),co3*co2*si1+si3*co1)
  ELSE IF ( j == 3 ) THEN
    result = equal(rmat(1,3),-co3*si2)
  END IF
ELSE IF ( i == 2 ) THEN
  IF ( j == 1 ) THEN
    result = equal(rmat(2,1),-si3*co2*co1-co3*si1)
  ELSE IF ( j == 2 ) THEN
    result = equal(rmat(2,2),-si3*co2*si1+co3*co1)
  ELSE IF ( j == 3 ) THEN
    result = equal(rmat(2,3),si3*si2)
  END IF
ELSE IF ( j == 1 ) THEN
  result = equal(rmat(3,1),si2*co1)
ELSE IF ( j == 2 ) THEN
  result = equal(rmat(3,2),si2*si1)
ELSE IF ( j == 3 ) THEN
  result = equal(rmat(3,3),co2)
END IF
checkrmat = result
END FUNCTION checkrmat
