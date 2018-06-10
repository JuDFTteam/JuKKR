subroutine rinvgj(ainv, a, arraydim, n)
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
  implicit none

! Dummy arguments
  integer :: arraydim, n
  real *8 :: a(arraydim, arraydim), ainv(arraydim, arraydim)

! Local variables
  integer :: icol, l, ll
  real *8 :: t, t1

  ainv(1, 1) = 0d0
!                                                        scan columns
  do icol = 1, n

!                                               make A(ICOL,ICOL) = 1
    t1 = 1.0d0/a(icol, icol)
    do l = (icol+1), n
      a(icol, l) = a(icol, l)*t1
    end do

    do l = 1, (icol-1)
      ainv(icol, l) = ainv(icol, l)*t1
    end do
    ainv(icol, icol) = t1

!                                    make A(LL,ICOL) = 0 for LL<>ICOL
    do ll = 1, n
      if (ll/=icol) then
        t = a(ll, icol)
        do l = (icol+1), n
          a(ll, l) = a(ll, l) - a(icol, l)*t
        end do

        do l = 1, (icol-1)
          ainv(ll, l) = ainv(ll, l) - ainv(icol, l)*t
        end do
        ainv(ll, icol) = -t1*t
      end if
    end do
  end do

end subroutine
