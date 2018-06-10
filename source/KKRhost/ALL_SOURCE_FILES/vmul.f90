! ************************************************************************
SUBROUTINE vmul(a,b,c)
! ************************************************************************

DOUBLE PRECISION, INTENT(IN)             :: a(*)
DOUBLE PRECISION, INTENT(IN)             :: b
DOUBLE PRECISION, INTENT(OUT)            :: c(*)

INTEGER :: i

DO  i=1,3
  c(i)=b*a(i)
END DO
RETURN
END SUBROUTINE vmul
