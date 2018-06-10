! ************************************************************************
SUBROUTINE veq(a,b)
! ************************************************************************
implicit none

DOUBLE PRECISION :: a(*)
DOUBLE PRECISION :: b(*)

INTEGER :: i

DO  i=1,3
  b(i)=a(i)
END DO
END SUBROUTINE veq
