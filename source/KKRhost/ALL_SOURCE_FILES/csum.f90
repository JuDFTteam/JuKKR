! 19.10.95 *************************************************************
COMPLEX*16 FUNCTION csum(n,v,iv)
! **********************************************************************
!        sum up the first N elements of the double complex
!        array V(*) with a stepwidth of IV
! ----------------------------------------------------------------------
!.. Scalar Arguments ..
      INTEGER IV,N
!..
!.. Array Arguments ..
      DOUBLE COMPLEX V(*)
!..
!.. Local Scalars ..
      DOUBLE COMPLEX VSUM
      INTEGER I,IBOT,ITOP
!..
IF (iv >= 0) THEN
  ibot = 1
  itop = 1 + (n-1)*iv
  
ELSE
  ibot = 1 - (n-1)*iv
  itop = 1
END IF

vsum = (0D0,0D0)
DO  i = ibot,itop,iv
  vsum = vsum + v(i)
END DO
csum = vsum
RETURN
END FUNCTION csum
