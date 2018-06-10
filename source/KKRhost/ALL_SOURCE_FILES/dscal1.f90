SUBROUTINE dscal1(n,da,dx,incx)
!- Scales a vector by a constant  dx(i) -> a * dx(i)
! ----------------------------------------------------------------------
!i Inputs:
!i   n     :lenght of dx and dy
!i   da    :constant
!i   dx    :vector
!i   incx  :incrementation for x
!o Outputs:
!o   dx    :vector
!r Remarks:
!r   Adapted from: jack dongarra, linpack, 3/11/78.
! ----------------------------------------------------------------------

implicit none 
! Passed parameters:                                                    
double precision da,dx(*) 
integer incx,n 
! Local parameters:                                                     
integer i,m,mp1,nincx 
!
IF( n <= 0 .OR. incx <= 0 )RETURN
IF(incx /= 1) THEN
! ----- code for increment not equal to 1
  nincx = n*incx
  DO i = 1,nincx,incx
    dx(i) = da*dx(i)
  END DO
ELSE
! ----- code for increment equal to 1
  m = MOD(n,5)
  IF( m /= 0 ) THEN
    DO i = 1,m
      dx(i) = da*dx(i)
    END DO
    IF( n < 5 ) RETURN
  END IF
  mp1 = m + 1
  DO i = mp1,n,5
    dx(i) = da*dx(i)
    dx(i + 1) = da*dx(i + 1)
    dx(i + 2) = da*dx(i + 2)
    dx(i + 3) = da*dx(i + 3)
    dx(i + 4) = da*dx(i + 4)
  END DO
END IF
END SUBROUTINE dscal1
