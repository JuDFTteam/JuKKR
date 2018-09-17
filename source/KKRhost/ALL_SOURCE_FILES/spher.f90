module mod_spher
  use :: mod_datatypes, only: dp
  private :: dp

contains

  subroutine spher(ylm, l, x)
    ! spherical harmonics except the facter exp(i*m*phi)

    ! m=-l to l , for given l.
    ! x=cos(theta)
    ! .. Scalar Arguments ..
    real (kind=dp) :: x
    integer :: l
    ! ..
    ! .. Array Arguments ..
    real (kind=dp) :: ylm(*)
    ! ..
    ! .. Local Scalars ..
    real (kind=dp) :: fac, ovr1, pi, qq
    integer :: i, ii, l2, lm, ln, m, nn
    ! ..
    ! .. Intrinsic Functions ..
    intrinsic :: abs, atan, real, sqrt
    ! ..
    pi = 4.0e0_dp*atan(1.0e0_dp)


    ovr1 = abs(x) - 1.e0_dp
    if (ovr1>0.1e-12_dp) then
      write (6, fmt=100) x
      stop
    else if (abs(ovr1)<1.e-10_dp) then
      if (x>0.0e0_dp) then
        fac = 1.0e0_dp
      else
        fac = (-1)**l
      end if
      l2 = 2*l + 1
      do i = 1, l2
        ylm(i) = 0.0e0_dp
      end do
      ylm(l+1) = sqrt(real(l2,kind=dp)/(4.0e0_dp*pi))*fac
      return
    end if

    ! l<0
    if (l<0) then
      write (6, fmt=*) ' === l=', l, ' < 0  : in sub.spher. ==='
      stop '=== stop in sub.spher. (l<0) ==='
      ! l=0
    else if (l==0) then
      ylm(1) = sqrt(1.0e0_dp/(4.0e0_dp*pi))
      ! l=1
    else if (l==1) then
      fac = sqrt(3.0e0_dp/(4.0e0_dp*pi))
      ylm(1) = fac*sqrt((1.0e0_dp-x*x)/2.0e0_dp)
      ylm(2) = fac*x
      ylm(3) = -ylm(1)
      ! l>1
    else
      ylm(1) = 1.0e0_dp
      ylm(2) = x
      do i = 2, l
        ylm(i+1) = ((2*i-1)*x*ylm(i)-(i-1)*ylm(i-1))/i
      end do
      fac = 1.0e0_dp/sqrt(1.0e0_dp-x*x)
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
      fac = sqrt((2*l+1)/(4.0e0_dp*pi))
      ylm(l+1) = fac*ylm(l+1)
      do m = 1, l
        fac = -fac/sqrt(real((l+m)*(l-m+1),kind=dp))
        lm = l + 1 + m
        ln = l + 1 - m
        qq = ylm(lm)
        ylm(lm) = fac*qq
        ylm(ln) = abs(fac)*qq
      end do
    end if

    return
100 format (/, /, 3x, '==invalid argument for spher; x=', d24.16, ' ==')
  end subroutine spher

end module mod_spher
