DOUBLE PRECISION FUNCTION ddot1(n,dx,incx,dy,incy)
!- Forms the dot product of two vectors.
! ----------------------------------------------------------------------
!i Inputs:
!i   n     :lenght of dx and dy
!i   dx    :first vector to mutiply
!i   incx  :incrementation for x
!i   dy    :second vector to mutiply
!i   incy  :incrementation for y
!o Outputs:
!o   ddot  :dot product of two vectors
!r Remarks:
!r    Adapted from: jack dongarra, linpack, 3/11/78.
! ----------------------------------------------------------------------
      implicit none
! Passed parameters:
      integer incx,incy,n
      double precision dx(*),dy(*)
! Local parameters:
      double precision dtemp
      integer i,ix,iy,m,mp1
!
ddot1 = 0.0D0
dtemp = 0.0D0
IF(n <= 0)RETURN
IF(incx /= 1.OR.incy /= 1)THEN
! ----- code for unequal increments or equal increments not equal to 1
  ix = 1
  iy = 1
  IF(incx < 0)ix = (-n+1)*incx + 1
  IF(incy < 0)iy = (-n+1)*incy + 1
  DO i = 1,n
    dtemp = dtemp + dx(ix)*dy(iy)
    ix = ix + incx
    iy = iy + incy
  END DO
  ddot1 = dtemp
ELSE
! ----- code for both increments equal to 1
  m = MOD(n,5)
  IF( m /= 0 ) THEN
    DO i = 1,m
      dtemp = dtemp + dx(i)*dy(i)
    END DO
    IF( n < 5 ) GO TO 60
  END IF
  mp1 = m + 1
  DO i = mp1,n,5
    dtemp = dtemp + dx(i)*dy(i) + dx(i+1)*dy(i+1) +                   &  &
        dx(i+2)*dy(i+2) + dx(i+3)*dy(i+3) + dx(i+4)*dy(i+4)
  END DO
  60   ddot1 = dtemp
END IF
END FUNCTION ddot1
