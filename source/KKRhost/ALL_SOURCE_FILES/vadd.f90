! ************************************************************************
SUBROUTINE vadd(a,b,c)
! ************************************************************************

DOUBLE PRECISION, INTENT(IN)             :: a(*)
DOUBLE PRECISION, INTENT(IN)             :: b(*)
DOUBLE PRECISION, INTENT(OUT)            :: c(*)

INTEGER :: i

DO  i=1,3
  c(i)=a(i)+b(i)
END DO
RETURN
END SUBROUTINE vadd
