SUBROUTINE simp3(f,fint,istart,iend,drdi)
!-----------------------------------------------------------------------
!     this subroutine does an integration from istart to iend of
!     the real function f with an extended 3-point-simpson :

!                          r(istart)

!                       fint = { f(r') dr'

!                           r(iend)

!-----------------------------------------------------------------------
!.. Scalar Arguments ..
DOUBLE PRECISION FINT
INTEGER IEND,ISTART
!..
!.. Array Arguments ..
DOUBLE PRECISION DRDI(*),F(*)
!..
!.. Local Scalars ..
DOUBLE PRECISION A1,A2
INTEGER I,IST
!..
!.. Intrinsic Functions ..
INTRINSIC MOD
!..
a1 = 4.0D0/3.0D0
a2 = 2.0D0/3.0D0

!---> initialize fint

IF (MOD(iend-istart,2) == 0) THEN
  fint = f(istart)*drdi(istart)/3.0D0
  ist = istart + 1
  
ELSE
  fint = (f(istart+3)*drdi(istart+3)- 5.0D0*f(istart+2)*drdi(istart+2)+  &
      19.0D0*f(istart+1)*drdi(istart+1)+ 9.0D0*f(istart)*drdi(istart))/24.0D0 +  &
      f(istart+1)*drdi(istart+1)/3.0D0
  ist = istart + 2
END IF

!---> calculate with an extended 3-point-simpson

DO  i = ist,iend - 1,2
  fint = fint + a1*f(i)*drdi(i) + a2*f(i+1)*drdi(i+1)
END DO
fint = fint - f(iend)*drdi(iend)/3.0D0

END SUBROUTINE simp3
