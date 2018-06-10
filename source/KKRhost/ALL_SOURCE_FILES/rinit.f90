SUBROUTINE rinit(n,a)
! **********************************************************************
! * Setting the first N values of a double precision array A to zero   *
! **********************************************************************
!..
!.. Arguments ..
      INTEGER N
      DOUBLE PRECISION A(*)
!..
!.. Locals ..
      INTEGER I,M,MP1
      DOUBLE PRECISION DZERO
!..
      DATA DZERO / 0.0D0 /
!..
!..
m = MOD(n,5)
IF ( m /= 0 ) THEN
  DO i = 1,m
    a(i) = dzero
  END DO
  IF ( n < 5 ) RETURN
END IF
mp1 = m + 1
DO i = mp1,n,5
  a(i  ) = dzero
  a(i+1) = dzero
  a(i+2) = dzero
  a(i+3) = dzero
  a(i+4) = dzero
END DO

END SUBROUTINE rinit
