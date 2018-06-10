subroutine idreals(darry, narry, iprint)
  implicit none

! PARAMETER definitions
  integer :: nsqr, nmul, divmax
  parameter (nsqr=7, nmul=5, divmax=15)
  double precision :: tol
  parameter (tol=1d-6)

! Dummy arguments
  integer :: iprint, narry
  double precision :: darry(narry)

! Local variables
  double precision :: dabs, dble, dsqrt, dsign
  integer :: div, i1, i2, idone(narry), imul(nmul), isqr(nsqr)
  double precision :: dsq, x, xn
  integer :: iabs, idnint

  data isqr/2, 3, 5, 6, 7, 8, 10/
  data imul/3, 7, 11, 13, 17/

! --> mark all numbers as unchecked

  do i1 = 1, narry
    idone(i1) = 0
  end do

! --> check darry**2/i integer?, i=1,divmax

  do div = 1, divmax
    dsq = dble(div)
    do i2 = 1, narry
      if (idone(i2)==0) then
        x = darry(i2)*darry(i2)*dsq
        xn = dnint(x)
        if (dabs(x-xn)/dsq<tol .and. xn/=0.d0) then
          if (iprint>4) write (1337, 100) dabs(darry(i2)), nint(x), div
          darry(i2) = dsign(1d0, darry(i2))*dsqrt(xn/dsq)
          idone(i2) = 1
        end if
      end if
    end do
  end do

! --> check darry/sqrt(n) =?=  i/j
!        n=2,3,5,6,7,8,10      i=1,divmax j=i*n

  do i1 = 1, nsqr
    do div = 1, divmax
      dsq = dsqrt(dble(div*div*isqr(i1)))
      do i2 = 1, narry
        if (idone(i2)==0) then
          x = darry(i2)*dsq
          xn = dnint(x)
          if (dabs(x-xn)/dsq<tol .and. xn/=0.d0) then
            if (iprint>4) write (1337, 110) dabs(darry(i2)), isqr(i1), &
              iabs(idnint(xn)), iabs(isqr(i1)*div)
            darry(i2) = xn/dsq
            idone(i2) = 1
          end if
        end if
      end do
    end do
  end do

! --> check darry = j/i * n ?
!        n=3,7,11,13,17

  do i1 = 1, nmul
    do div = 1, divmax
      dsq = dble(div*imul(i1))
      do i2 = 1, narry
        if (idone(i2)==0) then
          x = darry(i2)*dsq
          xn = dnint(x)
          if (dabs(x-xn)/dsq<tol .and. xn/=0.d0) then
            if (iprint>4) write (1337, 120) dabs(darry(i2)), imul(i1), &
              iabs(idnint(xn)), div
            darry(i2) = xn/dsq
            idone(i2) = 1
          end if
        end if
      end do
    end do
  end do
  return

100 format (8x, '< IDREALS > : identify ', f12.8, ' as dsqrt(', i3, '/', i3, &
    ')')
110 format (8x, '< IDREALS > : identify ', f12.8, ' as dsqrt(', i2, ')*', i3, &
    '/', i3)
120 format (8x, '< IDREALS > : identify ', f12.8, ' as 1/', i2, ' * ', i2, &
    '/', i1)

end subroutine
