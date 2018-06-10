SUBROUTINE spher(ylm,l,x)
!      spherical harmonics except the facter exp(i*m*phi)

!      m=-l to l , for given l.
!      x=cos(theta)
!     .. Scalar Arguments ..
      DOUBLE PRECISION X
      INTEGER L
!..
!.. Array Arguments ..
      DOUBLE PRECISION YLM(*)
!..
!.. Local Scalars ..
      DOUBLE PRECISION FAC,OVR1,PI,QQ
      INTEGER I,II,L2,LM,LN,M,NN
!..
!.. Intrinsic Functions ..
      INTRINSIC ABS,ATAN,DBLE,SQRT
!     ..
pi = 4.0D0*ATAN(1.0D0)


ovr1 = ABS(x) - 1.d0
IF (ovr1 > 0.1D-12) THEN
  WRITE (6,FMT=9000) x
  STOP
ELSE IF (ABS(ovr1) < 1.d-10) THEN
  IF (x > 0.0D0) THEN
    fac = 1.0D0
  ELSE
    fac = (-1)**l
  END IF
  l2 = 2*l + 1
  DO  i = 1,l2
    ylm(i) = 0.0D0
  END DO
  ylm(l+1) = SQRT(DBLE(l2)/ (4.0D0*pi))*fac
  RETURN
END IF

! l<0
IF (l < 0) THEN
  WRITE (6,FMT=*) ' === l=',l,' < 0  : in sub.spher. ==='
  STOP '=== stop in sub.spher. (l<0) ==='
! l=0
ELSE IF (l == 0) THEN
  ylm(1) = SQRT(1.0D0/ (4.0D0*pi))
! l=1
ELSE IF (l == 1) THEN
  fac = SQRT(3.0D0/ (4.0D0*pi))
  ylm(1) = fac*SQRT((1.0D0-x*x)/2.0D0)
  ylm(2) = fac*x
  ylm(3) = -ylm(1)
! l>1
ELSE
  ylm(1) = 1.0D0
  ylm(2) = x
  DO  i = 2,l
    ylm(i+1) = ((2*i-1)*x*ylm(i)- (i-1)*ylm(i-1))/i
  END DO
  fac = 1.0D0/SQRT(1.0D0-x*x)
  DO  m = 1,l
    lm = l + m
    ylm(lm+1) = fac* (- (l-m+1)*x*ylm(lm)+ (lm-1)*ylm(l))
    IF (m < l) THEN
      nn = m + 1
      DO  i = nn,l
        ii = l - i + nn
        ylm(ii) = fac* (- (ii-m)*x*ylm(ii)+ (ii+m-2)*ylm(ii-1))
      END DO
    END IF
  END DO
  fac = SQRT((2*l+1)/ (4.0D0*pi))
  ylm(l+1) = fac*ylm(l+1)
  DO  m = 1,l
    fac = -fac/SQRT(DBLE((l+m)* (l-m+1)))
    lm = l + 1 + m
    ln = l + 1 - m
    qq = ylm(lm)
    ylm(lm) = fac*qq
    ylm(ln) = ABS(fac)*qq
  END DO
END IF

RETURN
9000 FORMAT (/,/,3X,'==invalid argument for spher; x=',d24.16,' ==')
END SUBROUTINE spher
