FUNCTION ylag(xi,x,y,ind1,n1,imax)
!   ********************************************************************
!   *                                                                  *
!   * lagrangian interpolation                                         *
!   * xi is interpolated entry into x-array                            *
!   * n is the order of lagrangran interpolation                       *
!   * y is array from which ylag is obtained by interpolation          *
!   * ind is the min-i for x(i).gt.xi                                  *
!   * if ind=0,x-array will be searched                                *
!   * imax is max index of x-and y-arrays                              *
!   *                                                                  *
!   * 07/12/94  HE  arg. IEX removed                                   *
!   ********************************************************************

IMPLICIT NONE

! Dummy arguments
INTEGER IMAX,IND1,N1
REAL*8 XI
REAL*8 X(IMAX),Y(IMAX)
REAL*8 YLAG

! Local variables
REAL*8 D,P,S,XD
INTEGER I,IND,INL,INU,J,N
SAVE D,I,IND,INL,INU,J,N,P,S,XD

ind = ind1
n = n1
IF ( n > imax ) n = imax
IF ( ind > 0 ) GO TO 200
DO j = 1,imax
  IF ( ABS(xi-x(j)) < 1.0D-12 ) GO TO 600
  IF ( xi < x(j) ) GO TO 100
  IF ( xi == x(j) ) GO TO 600
END DO
GO TO 300
100  CONTINUE
ind = j
200  CONTINUE
IF ( ind > 1 ) THEN
END IF
inl = ind - (n+1)/2
IF ( inl <= 0 ) inl = 1
inu = inl + n - 1
IF ( inu <= imax ) GO TO 400
300  CONTINUE
inl = imax - n + 1
inu = imax
400  CONTINUE
s = 0.0D0
p = 1.0D0
DO j = inl,inu
  p = p*(xi-x(j))
  d = 1.0D0
  DO i = inl,inu
    IF ( i /= j ) THEN
      xd = x(j)
    ELSE
      xd = xi
    END IF
    d = d*(xd-x(i))
  END DO
  s = s + y(j)/d
END DO
ylag = s*p
500  CONTINUE
RETURN
600  CONTINUE
ylag = y(j)
GO TO 500
END FUNCTION ylag


