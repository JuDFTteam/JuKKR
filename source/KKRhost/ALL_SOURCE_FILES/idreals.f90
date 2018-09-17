module mod_idreals
  use :: mod_datatypes, only: dp
  private :: dp

contains

  ! < check iif entries of darry are fractions or roots of common numbers and
  ! replace accordingly for hogher accuracy
  subroutine idreals(darry, narry, iprint)
    implicit none

    ! PARAMETER definitions
    integer, parameter :: nsqr = 7
    integer, parameter :: nmul = 5
    integer, parameter :: divmax = 15
    real (kind=dp), parameter :: tol = 1e-6_dp
    real (kind=dp), parameter :: eps = 1e-14_dp

    ! Dummy arguments
    integer :: iprint, narry
    real (kind=dp) :: darry(narry)

    ! Local variables
    integer :: div, i1, i2, idone(narry), imul(nmul), isqr(nsqr)
    real (kind=dp) :: dsq, x, xn

    data isqr/2, 3, 5, 6, 7, 8, 10/
    data imul/3, 7, 11, 13, 17/

    ! --> mark all numbers as unchecked

    do i1 = 1, narry
      idone(i1) = 0
    end do

    ! --> check darry**2/i integer?, i=1,divmax

    do div = 1, divmax
      dsq = real(div, kind=dp)
      do i2 = 1, narry
        if (idone(i2)==0) then
          x = darry(i2)*darry(i2)*dsq
          xn = nint(x)
          if (abs(x-xn)/dsq<tol .and. abs(xn)>eps) then
            if (iprint>4) write (1337, 100) abs(darry(i2)), nint(x), div
            darry(i2) = sign(1e0_dp, darry(i2))*sqrt(xn/dsq)
            idone(i2) = 1
          end if
        end if
      end do
    end do

    ! --> check darry/sqrt(n) =?=  i/j
    ! n=2,3,5,6,7,8,10      i=1,divmax j=i*n

    do i1 = 1, nsqr
      do div = 1, divmax
        dsq = sqrt(real(div*div*isqr(i1),kind=dp))
        do i2 = 1, narry
          if (idone(i2)==0) then
            x = darry(i2)*dsq
            xn = nint(x)
            if (abs(x-xn)/dsq<tol .and. abs(xn)>eps) then
              if (iprint>4) write (1337, 110) abs(darry(i2)), isqr(i1), abs(nint(xn)), abs(isqr(i1)*div)
              darry(i2) = xn/dsq
              idone(i2) = 1
            end if
          end if
        end do
      end do
    end do

    ! --> check darry = j/i * n ?
    ! n=3,7,11,13,17

    do i1 = 1, nmul
      do div = 1, divmax
        dsq = real(div*imul(i1), kind=dp)
        do i2 = 1, narry
          if (idone(i2)==0) then
            x = darry(i2)*dsq
            xn = nint(x)
            if (abs(x-xn)/dsq<tol .and. abs(xn)>eps) then
              if (iprint>4) write (1337, 120) abs(darry(i2)), imul(i1), abs(nint(xn)), div
              darry(i2) = xn/dsq
              idone(i2) = 1
            end if
          end if
        end do
      end do
    end do
    return

100 format (8x, '< IDREALS > : identify ', f12.8, ' as sqrt(', i3, '/', i3, ')')
110 format (8x, '< IDREALS > : identify ', f12.8, ' as sqrt(', i2, ')*', i3, '/', i3)
120 format (8x, '< IDREALS > : identify ', f12.8, ' as 1/', i2, ' * ', i2, '/', i1)

  end subroutine idreals

end module mod_idreals
