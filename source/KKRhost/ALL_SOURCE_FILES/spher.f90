subroutine spher(ylm, l, x)
!      spherical harmonics except the facter exp(i*m*phi)

!      m=-l to l , for given l.
!      x=cos(theta)
!     .. Scalar Arguments ..
  double precision :: x
  integer :: l
!..
!.. Array Arguments ..
  double precision :: ylm(*)
!..
!.. Local Scalars ..
  double precision :: fac, ovr1, pi, qq
  integer :: i, ii, l2, lm, ln, m, nn
!..
!.. Intrinsic Functions ..
  intrinsic :: abs, atan, dble, sqrt
!     ..
  pi = 4.0d0*atan(1.0d0)


  ovr1 = abs(x) - 1.d0
  if (ovr1>0.1d-12) then
    write (6, fmt=100) x
    stop
  else if (abs(ovr1)<1.d-10) then
    if (x>0.0d0) then
      fac = 1.0d0
    else
      fac = (-1)**l
    end if
    l2 = 2*l + 1
    do i = 1, l2
      ylm(i) = 0.0d0
    end do
    ylm(l+1) = sqrt(dble(l2)/(4.0d0*pi))*fac
    return
  end if

! l<0
  if (l<0) then
    write (6, fmt=*) ' === l=', l, ' < 0  : in sub.spher. ==='
    stop '=== stop in sub.spher. (l<0) ==='
! l=0
  else if (l==0) then
    ylm(1) = sqrt(1.0d0/(4.0d0*pi))
! l=1
  else if (l==1) then
    fac = sqrt(3.0d0/(4.0d0*pi))
    ylm(1) = fac*sqrt((1.0d0-x*x)/2.0d0)
    ylm(2) = fac*x
    ylm(3) = -ylm(1)
! l>1
  else
    ylm(1) = 1.0d0
    ylm(2) = x
    do i = 2, l
      ylm(i+1) = ((2*i-1)*x*ylm(i)-(i-1)*ylm(i-1))/i
    end do
    fac = 1.0d0/sqrt(1.0d0-x*x)
    do m = 1, l
      lm = l + m
      ylm(lm+1) = fac*(-(l-m+1)*x*ylm(lm)+(lm-1)*ylm(l))
      if (m<l) then
        nn = m + 1
        do i = nn, l
          ii = l - i + nn
          ylm(ii) = fac*(-(ii-m)*x*ylm(ii)+(ii+m-2)*ylm(ii-1))
        end do
      end if
    end do
    fac = sqrt((2*l+1)/(4.0d0*pi))
    ylm(l+1) = fac*ylm(l+1)
    do m = 1, l
      fac = -fac/sqrt(dble((l+m)*(l-m+1)))
      lm = l + 1 + m
      ln = l + 1 - m
      qq = ylm(lm)
      ylm(lm) = fac*qq
      ylm(ln) = abs(fac)*qq
    end do
  end if

  return
100 format (/, /, 3x, '==invalid argument for spher; x=', d24.16, ' ==')
end subroutine
