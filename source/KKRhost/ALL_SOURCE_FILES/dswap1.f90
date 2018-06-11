    Subroutine dswap1(n, dx, incx, dy, incy)
      Use mod_datatypes, Only: dp
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

      Implicit None
! Passed parameters:                                                    
      Integer :: incx, incy, n
      Real (Kind=dp) :: dx(*), dy(*)
! Local parameters:                                                     
      Real (Kind=dp) :: dtemp
      Integer :: i, ix, iy, m, mp1

      If (n<=0) Return
      If (incx/=1 .Or. incy/=1) Then
! ----- code for unequal increments or equal increments not equal to 1
        ix = 1
        iy = 1
        If (incx<0) ix = (-n+1)*incx + 1
        If (incy<0) iy = (-n+1)*incy + 1
        Do i = 1, n
          dtemp = dx(ix)
          dx(ix) = dy(iy)
          dy(iy) = dtemp
          ix = ix + incx
          iy = iy + incy
        End Do
      Else
! ----- code for both increments equal to 1
        m = mod(n, 3)
        If (m/=0) Then
          Do i = 1, m
            dtemp = dx(i)
            dx(i) = dy(i)
            dy(i) = dtemp
          End Do
          If (n<3) Return
        End If
        mp1 = m + 1
        Do i = mp1, n, 3
          dtemp = dx(i)
          dx(i) = dy(i)
          dy(i) = dtemp
          dtemp = dx(i+1)
          dx(i+1) = dy(i+1)
          dy(i+1) = dtemp
          dtemp = dx(i+2)
          dx(i+2) = dy(i+2)
          dy(i+2) = dtemp
        End Do
      End If
    End Subroutine
