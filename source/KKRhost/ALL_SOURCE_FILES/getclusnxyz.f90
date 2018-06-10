SUBROUTINE getclusnxyz(clurad,bravais,ndim,cluradsq,nbr)
! **********************************************************************
! *                                                                    *
! * Given a spherical cluster of radius CLURAD it determines the three *
! * integers N1,N2,N3 such that any vector                             *
! *                                                                    *
! *    R_i = r_i + SUM_j  N_j * a_j                                    *
! *                                                                    *
! *  with i = 1,NAEZ and a_j the primitive Bravais vectors, is inside  *
! *  the cluster. Subroutine also returns the CLURAD**2 value          *
! *                                                                    *
! **********************************************************************
      IMPLICIT NONE
! ..
! ..  Arguments
      DOUBLE PRECISION CLURAD,CLURADSQ
      DOUBLE PRECISION BRAVAIS(3,3)
      INTEGER NDIM,NBR(3)
! .. 
! ..  Locals
      DOUBLE PRECISION DR(3)
      INTEGER I,J
      INTEGER INT
! ..
! ..
DO i = 1,ndim
  dr(i) = 0D0
  DO j = 1,ndim
    dr(i) = dr(i) + bravais(j,i)*bravais(j,i)
  END DO
  dr(i) = SQRT(dr(i))
END DO

IF ( ABS(clurad) < 1D-6 ) THEN
  DO i = 1,ndim
    nbr(i) = 0
  END DO
  cluradsq = 1D10
ELSE
  DO i = 1,ndim
    nbr(i) = INT(clurad/dr(i)) + 2
  END DO
  cluradsq = clurad*clurad
END IF

END SUBROUTINE getclusnxyz
