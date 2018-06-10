SUBROUTINE idreals(darry,narry,iprint)
IMPLICIT NONE

! PARAMETER definitions
INTEGER NSQR,NMUL,DIVMAX
PARAMETER (NSQR=7,NMUL=5,DIVMAX=15)
DOUBLE PRECISION TOL
PARAMETER (TOL=1D-6)

! Dummy arguments
INTEGER IPRINT,NARRY
DOUBLE PRECISION DARRY(NARRY)

! Local variables
DOUBLE PRECISION DABS,DBLE,DSQRT,DSIGN
INTEGER DIV,I1,I2,IDONE(NARRY),IMUL(NMUL),ISQR(NSQR)
DOUBLE PRECISION DSQ,X,XN
INTEGER IABS,IDNINT

DATA ISQR/2,3,5,6,7,8,10/
DATA IMUL/3,7,11,13,17/

! --> mark all numbers as unchecked

DO i1 = 1,narry
  idone(i1) = 0
END DO

! --> check darry**2/i integer?, i=1,divmax

DO div = 1,divmax
  dsq = DBLE(div)
  DO i2 = 1,narry
    IF ( idone(i2) == 0 ) THEN
      x = darry(i2)*darry(i2)*dsq
      xn = DNINT(x)
      IF ( DABS(x-xn)/dsq < tol .AND. xn /= 0.d0 ) THEN
        IF (iprint > 4) WRITE (1337,99000) DABS(darry(i2)),nint(x),div
        darry(i2) = DSIGN(1D0,darry(i2))*DSQRT(xn/dsq)
        idone(i2) = 1
      END IF
    END IF
  END DO
END DO

! --> check darry/sqrt(n) =?=  i/j
!        n=2,3,5,6,7,8,10      i=1,divmax j=i*n

DO i1 = 1,nsqr
  DO div = 1,divmax
    dsq = DSQRT(DBLE(div*div*isqr(i1)))
    DO i2 = 1,narry
      IF ( idone(i2) == 0 ) THEN
        x = darry(i2)*dsq
        xn = DNINT(x)
        IF ( DABS(x-xn)/dsq < tol .AND. xn /= 0.d0 ) THEN
          IF (iprint > 4) WRITE (1337,99001) DABS(darry(i2)),isqr(i1),  &
              IABS(IDNINT(xn)),IABS(isqr(i1)*div)
          darry(i2) = xn/dsq
          idone(i2) = 1
        END IF
      END IF
    END DO
  END DO
END DO

! --> check darry = j/i * n ?
!        n=3,7,11,13,17

DO i1 = 1,nmul
  DO div = 1,divmax
    dsq = DBLE(div*imul(i1))
    DO i2 = 1,narry
      IF ( idone(i2) == 0 ) THEN
        x = darry(i2)*dsq
        xn = DNINT(x)
        IF ( DABS(x-xn)/dsq < tol .AND. xn /= 0.d0 ) THEN
          IF (iprint > 4) WRITE(1337,99002)  &
              DABS(darry(i2)),imul(i1),IABS(IDNINT(xn)),div
          darry(i2) = xn/dsq
          idone(i2) = 1
        END IF
      END IF
    END DO
  END DO
END DO
RETURN

99000 FORMAT (8X,'< IDREALS > : identify ',f12.8, ' as dsqrt(',i3,'/',i3,')')
99001 FORMAT (8X,'< IDREALS > : identify ',f12.8,  &
    ' as dsqrt(',i2,')*',i3,'/',i3)
99002 FORMAT (8X,'< IDREALS > : identify ',f12.8,' as 1/',i2, ' * ',i2,'/',i1)

END SUBROUTINE idreals
