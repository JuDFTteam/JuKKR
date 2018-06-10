DOUBLE PRECISION FUNCTION ssum(n,v,iv)
! **********************************************************************
!        sum up the first N elements of the double precision
!        array V(*) with a stepwidth of IV
! ----------------------------------------------------------------------
implicit none
!.. Scalar Arguments ..
INTEGER IV,N
!..
!.. Array Arguments ..
DOUBLE PRECISION V(*)
!..
!.. Local Scalars ..
DOUBLE PRECISION VSUM
INTEGER I,IBOT,ITOP
!..
IF (iv >= 0) THEN
  ibot = 1
  itop = 1 + (n-1)*iv
  
ELSE
  ibot = 1 - (n-1)*iv
  itop = 1
END IF

vsum = 0.0D0
DO  i = ibot,itop,iv
  vsum = vsum + v(i)
END DO
ssum = vsum
END FUNCTION ssum
