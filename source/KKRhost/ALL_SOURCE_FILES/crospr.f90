! ************************************************************************
SUBROUTINE crospr(x,y,z)
! ************************************************************************
!     CROSP COMPUTES THE CROSS PRODUCT OF X AND Y RETURNING
!     IT INTO Z.
! ------------------------------------------------------------------------
implicit none
DOUBLE PRECISION, INTENT(IN)             :: x(*)
DOUBLE PRECISION, INTENT(IN)             :: y(*)
DOUBLE PRECISION, INTENT(OUT)            :: z(*)

z(1)=x(2)*y(3)-x(3)*y(2)
z(2)=x(3)*y(1)-x(1)*y(3)
z(3)=x(1)*y(2)-x(2)*y(1)
RETURN
END SUBROUTINE crospr
