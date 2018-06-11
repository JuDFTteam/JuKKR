    Function ddot1(n, dx, incx, dy, incy)
      Use mod_datatypes, Only: dp
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
      Implicit None
      Real (Kind=dp) :: ddot1
! Passed parameters:
      Integer :: incx, incy, n
      Real (Kind=dp) :: dx(*), dy(*)
! Local parameters:
      Real (Kind=dp) :: dtemp
      Integer :: i, ix, iy, m, mp1
!
      ddot1 = 0.0E0_dp
      dtemp = 0.0E0_dp
      If (n<=0) Return
      If (incx/=1 .Or. incy/=1) Then
! ----- code for unequal increments or equal increments not equal to 1
        ix = 1
        iy = 1
        If (incx<0) ix = (-n+1)*incx + 1
        If (incy<0) iy = (-n+1)*incy + 1
        Do i = 1, n
          dtemp = dtemp + dx(ix)*dy(iy)
          ix = ix + incx
          iy = iy + incy
        End Do
        ddot1 = dtemp
      Else
! ----- code for both increments equal to 1
        m = mod(n, 5)
        If (m/=0) Then
          Do i = 1, m
            dtemp = dtemp + dx(i)*dy(i)
          End Do
          If (n<5) Go To 100
        End If
        mp1 = m + 1
        Do i = mp1, n, 5
          dtemp = dtemp + dx(i)*dy(i) + dx(i+1)*dy(i+1) + dx(i+2)*dy(i+2) + &
            dx(i+3)*dy(i+3) + dx(i+4)*dy(i+4)
        End Do
100     ddot1 = dtemp
      End If
    End Function
