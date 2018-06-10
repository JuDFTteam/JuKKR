SUBROUTINE cmatmul(n,m,a,b,c)
!   ********************************************************************
!   *                                                                  *
!   *   perform  the matrix-matrix operation           C = A * B       *
!   *                                                                  *
!   *   A,B,C   complex  SQUARE  N x N - matrices                      *
!   *   N       dimension of A, B and C                                *
!   *   M       array size of A, B, C with M >= N                      *
!   *                                                                  *
!   ********************************************************************

IMPLICIT none

! PARAMETER definitions
DOUBLE COMPLEX C0
PARAMETER (C0=(0.0D0,0.0D0))

! Dummy arguments
INTEGER M,N
DOUBLE COMPLEX A(M,M),B(M,M),C(M,M)

! Local variables
DOUBLE COMPLEX BLJ
INTEGER I,J,L

DO j = 1,n
  DO i = 1,n
    c(i,j) = c0
  END DO
END DO

DO j = 1,n
  DO l = 1,n
    blj = b(l,j)
    IF ( blj /= c0 ) THEN
      DO i = 1,n
        c(i,j) = c(i,j) + a(i,l)*blj
      END DO
    END IF
  END DO
END DO

END SUBROUTINE cmatmul
