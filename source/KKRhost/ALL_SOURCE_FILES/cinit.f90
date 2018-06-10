SUBROUTINE cinit(n,a)
! **********************************************************************
! * Setting the first N values of a double complex array A to zero     *
! **********************************************************************
!     ..
!     .. Arguments ..
      INTEGER N
      DOUBLE COMPLEX A(*)
!     ..
!     .. Locals
      INTEGER I,M,MP1
      DOUBLE COMPLEX CZERO
!     ..
      DATA CZERO / (0.0D0,0.0D0) /
!     ..
m = MOD(n,5)
IF ( m /= 0 ) THEN
  DO i = 1,m
    a(i) = czero
  END DO
  IF ( n < 5 ) RETURN
END IF
mp1 = m + 1
DO i = mp1,n,5
  a(i  ) = czero
  a(i+1) = czero
  a(i+2) = czero
  a(i+3) = czero
  a(i+4) = czero
END DO

END SUBROUTINE cinit
