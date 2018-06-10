! ************************************************************************
SUBROUTINE scalpr(x,y,z)
! ************************************************************************
!     SCALSP COMPUTES THE scalar PRODUCT OF X AND Y RETURNING
!     IT INTO Z.
! ------------------------------------------------------------------------


DOUBLE PRECISION, INTENT(IN)             :: x(*)
DOUBLE PRECISION, INTENT(IN)             :: y(*)
DOUBLE PRECISION, INTENT(OUT)            :: z

z= x(1)*y(1)+x(2)*y(2)+x(3)*y(3)
RETURN
END SUBROUTINE scalpr
