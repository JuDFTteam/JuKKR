    Subroutine rinvgj(ainv, a, arraydim, n)
      Use mod_datatypes, Only: dp
!   ********************************************************************
!   *                                                                  *
!   *                      AINV = A**(-1)                              *
!   *                                                                  *
!   *  invert A using the GAUSS-JORDAN - algorithm                     *
!   *  the 1- matrix is not set up and use is made of its structure    *
!   *                                                                  *
!   *                    REAL*8 VERSION                                *
!   *                                                                  *
!   ********************************************************************
      Implicit None

! Dummy arguments
      Integer :: arraydim, n
      Real (Kind=dp) :: a(arraydim, arraydim), ainv(arraydim, arraydim)

! Local variables
      Integer :: icol, l, ll
      Real (Kind=dp) :: t, t1

      ainv(1, 1) = 0E0_dp
!                                                        scan columns
      Do icol = 1, n

!                                               make A(ICOL,ICOL) = 1
        t1 = 1.0E0_dp/a(icol, icol)
        Do l = (icol+1), n
          a(icol, l) = a(icol, l)*t1
        End Do

        Do l = 1, (icol-1)
          ainv(icol, l) = ainv(icol, l)*t1
        End Do
        ainv(icol, icol) = t1

!                                    make A(LL,ICOL) = 0 for LL<>ICOL
        Do ll = 1, n
          If (ll/=icol) Then
            t = a(ll, icol)
            Do l = (icol+1), n
              a(ll, l) = a(ll, l) - a(icol, l)*t
            End Do

            Do l = 1, (icol-1)
              ainv(ll, l) = ainv(ll, l) - ainv(icol, l)*t
            End Do
            ainv(ll, icol) = -t1*t
          End If
        End Do
      End Do

    End Subroutine
