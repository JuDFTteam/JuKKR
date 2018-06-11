    Subroutine dscal1(n, da, dx, incx)
      Use mod_datatypes, Only: dp
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

      Implicit None
! Passed parameters:                                                    
      Real (Kind=dp) :: da, dx(*)
      Integer :: incx, n
! Local parameters:                                                     
      Integer :: i, m, mp1, nincx
!
      If (n<=0 .Or. incx<=0) Return
      If (incx/=1) Then
! ----- code for increment not equal to 1
        nincx = n*incx
        Do i = 1, nincx, incx
          dx(i) = da*dx(i)
        End Do
      Else
! ----- code for increment equal to 1
        m = mod(n, 5)
        If (m/=0) Then
          Do i = 1, m
            dx(i) = da*dx(i)
          End Do
          If (n<5) Return
        End If
        mp1 = m + 1
        Do i = mp1, n, 5
          dx(i) = da*dx(i)
          dx(i+1) = da*dx(i+1)
          dx(i+2) = da*dx(i+2)
          dx(i+3) = da*dx(i+3)
          dx(i+4) = da*dx(i+4)
        End Do
      End If
    End Subroutine
