      SUBROUTINE RINVGJ(AINV,A,ADIM,N)
C   ********************************************************************
C   *                                                                  *
C   *                      AINV = A**(-1)                              *
C   *                                                                  *
C   *  invert A using the GAUSS-JORDAN - algorithm                     *
C   *  the 1- matrix is not set up and use is made of its structure    *
C   *                                                                  *
C   *                    double precision VERSION                      *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C
C Dummy arguments
C
      INTEGER, intent(in) :: ADIM, N
      double precision, intent(inout) ::  A(ADIM,ADIM)
      double precision, intent(out) :: AINV(ADIM,ADIM)
C
C Local variables
C
      INTEGER :: ICOL, L, LL
      double precision :: T, T1
C
      AINV(1,1) = 0.D0
C                                                        scan columns
      DO ICOL = 1,N
C
C                                               make A(ICOL,ICOL) = 1
         T1 = 1.0D0/A(ICOL,ICOL)
         DO L = (ICOL+1),N
            A(ICOL,L) = A(ICOL,L)*T1
         END DO
C
         DO L = 1,(ICOL-1)
            AINV(ICOL,L) = AINV(ICOL,L)*T1
         END DO
         AINV(ICOL,ICOL) = T1
C
C                                    make A(LL,ICOL) = 0 for LL<>ICOL
         DO LL = 1,N
            IF (LL /= ICOL) THEN
               T = A(LL,ICOL)
               DO L = (ICOL+1),N
                  A(LL,L) = A(LL,L) - A(ICOL,L)*T
               END DO
C
               DO L = 1,(ICOL-1)
                  AINV(LL,L) = AINV(LL,L) - AINV(ICOL,L)*T
               END DO
               AINV(LL,ICOL) = -T1*T
            END IF
         END DO
      END DO
C
      END
