SUBROUTINE dirbslag(xi,y1i,y2i,y3i,y4i,y1,y2,y3,y4,ind1,n,imax)
!   ********************************************************************
!   *                                                                  *
!   *      lagrangian interpolation of Y(X) at position XI             *
!   *                                                                  *
!   *      XI      entry into x-array                                  *
!   *              for regular solution:   X(IND1-1) < XI <=X(IND1)    *
!   *              for irregular solution: X(IND1)   < XI <=X(IND1+1)  *
!   *      X/Y     X/Y-arrays                                          *
!   *      N       order of lagrangian interpolation                   *
!   *      IND     min-I for which  X(I) > XI                          *
!   *      IMAX    max index of X/Y-arrays                             *
!   *                                                                  *
!   ********************************************************************
IMPLICIT NONE

! Dummy arguments
INTEGER IMAX,IND1,N
REAL*8 XI
REAL*8 Y1I,Y2I,Y3I,Y4I
REAL*8 Y1(IMAX),Y2(IMAX),Y3(IMAX),Y4(IMAX)

! Local variables
REAL*8 D,P,XD
INTEGER I,IND,INL,INU,J

ind = ind1
IF ( ABS(xi-DBLE(ind)) < 1.0D-12 ) THEN
  y1i = y1(ind)
  y2i = y2(ind)
  y3i = y3(ind)
  y4i = y4(ind)
  RETURN
END IF
! ------------------------------------- shift IND for irregular solution
IF ( xi > DBLE(ind) ) ind = ind + 1

inl = MAX(1,ind-(n+1)/2)
inu = inl + n - 1

IF ( inu > imax ) THEN
  inl = imax - n + 1
  inu = imax
END IF

y1i = 0.0D0
y2i = 0.0D0
y3i = 0.0D0
y4i = 0.0D0
p = 1.0D0
DO j = inl,inu
  p = p*(xi-DBLE(j))
  d = 1.0D0
  DO i = inl,inu
    IF ( i /= j ) THEN
      xd = DBLE(j)
    ELSE
      xd = xi
    END IF
    d = d*(xd-DBLE(i))
  END DO
  
  y1i = y1i + y1(j)/d
  y2i = y2i + y2(j)/d
  y3i = y3i + y3(j)/d
  y4i = y4i + y4(j)/d
  
END DO
y1i = y1i*p
y2i = y2i*p
y3i = y3i*p
y4i = y4i*p
END SUBROUTINE dirbslag
