SUBROUTINE dswap1 (n,dx,incx,dy,incy)
!-Interchanges two vectors
! ----------------------------------------------------------------------
!i Inputs:
!i   n     :lenght of dx and dy
!io  dx    :vector
!i   incx  :incrementation for x
!io  dy    :vector
!i   incy  :incrementation for y
!o Outputs:
!io  dx    :vector
!io  dy    :vector
!r Remarks:
!r Adapted from:  jack dongarra, linpack, 3/11/78.
! ----------------------------------------------------------------------

      implicit none 
! Passed parameters:                                                    
      integer incx,incy,n 
      double precision dx(*),dy(*) 
! Local parameters:                                                     
      double precision dtemp 
      integer i,ix,iy,m,mp1 
                                                                        
IF(n <= 0)RETURN
IF(incx /= 1.OR.incy /= 1) THEN
! ----- code for unequal increments or equal increments not equal to 1
  ix = 1
  iy = 1
  IF(incx < 0)ix = (-n+1)*incx + 1
  IF(incy < 0)iy = (-n+1)*incy + 1
  DO i = 1,n
    dtemp = dx(ix)
    dx(ix) = dy(iy)
    dy(iy) = dtemp
    ix = ix + incx
    iy = iy + incy
  END DO
ELSE
! ----- code for both increments equal to 1
  m = MOD(n,3)
  IF( m /= 0 )  THEN
    DO i = 1,m
      dtemp = dx(i)
      dx(i) = dy(i)
      dy(i) = dtemp
    END DO
    IF( n < 3 ) RETURN
  END IF
  mp1 = m + 1
  DO i = mp1,n,3
    dtemp = dx(i)
    dx(i) = dy(i)
    dy(i) = dtemp
    dtemp = dx(i + 1)
    dx(i + 1) = dy(i + 1)
    dy(i + 1) = dtemp
    dtemp = dx(i + 2)
    dx(i + 2) = dy(i + 2)
    dy(i + 2) = dtemp
  END DO
END IF
END SUBROUTINE dswap1
