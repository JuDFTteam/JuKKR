subroutine dswap1(n, dx, incx, dy, incy)
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
  integer :: incx, incy, n
  double precision :: dx(*), dy(*)
! Local parameters:                                                     
  double precision :: dtemp
  integer :: i, ix, iy, m, mp1

  if (n<=0) return
  if (incx/=1 .or. incy/=1) then
! ----- code for unequal increments or equal increments not equal to 1
    ix = 1
    iy = 1
    if (incx<0) ix = (-n+1)*incx + 1
    if (incy<0) iy = (-n+1)*incy + 1
    do i = 1, n
      dtemp = dx(ix)
      dx(ix) = dy(iy)
      dy(iy) = dtemp
      ix = ix + incx
      iy = iy + incy
    end do
  else
! ----- code for both increments equal to 1
    m = mod(n, 3)
    if (m/=0) then
      do i = 1, m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
      end do
      if (n<3) return
    end if
    mp1 = m + 1
    do i = mp1, n, 3
      dtemp = dx(i)
      dx(i) = dy(i)
      dy(i) = dtemp
      dtemp = dx(i+1)
      dx(i+1) = dy(i+1)
      dy(i+1) = dtemp
      dtemp = dx(i+2)
      dx(i+2) = dy(i+2)
      dy(i+2) = dtemp
    end do
  end if
end subroutine
