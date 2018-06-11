    Subroutine spher(ylm, l, x)
      Use mod_datatypes, Only: dp
!      spherical harmonics except the facter exp(i*m*phi)

!      m=-l to l , for given l.
!      x=cos(theta)
!     .. Scalar Arguments ..
      Real (Kind=dp) :: x
      Integer :: l
!..
!.. Array Arguments ..
      Real (Kind=dp) :: ylm(*)
!..
!.. Local Scalars ..
      Real (Kind=dp) :: fac, ovr1, pi, qq
      Integer :: i, ii, l2, lm, ln, m, nn
!..
!.. Intrinsic Functions ..
      Intrinsic :: abs, atan, real, sqrt
!     ..
      pi = 4.0E0_dp*atan(1.0E0_dp)


      ovr1 = abs(x) - 1.E0_dp
      If (ovr1>0.1E-12_dp) Then
        Write (6, Fmt=100) x
        Stop
      Else If (abs(ovr1)<1.E-10_dp) Then
        If (x>0.0E0_dp) Then
          fac = 1.0E0_dp
        Else
          fac = (-1)**l
        End If
        l2 = 2*l + 1
        Do i = 1, l2
          ylm(i) = 0.0E0_dp
        End Do
        ylm(l+1) = sqrt(real(l2,kind=dp)/(4.0E0_dp*pi))*fac
        Return
      End If

! l<0
      If (l<0) Then
        Write (6, Fmt=*) ' === l=', l, ' < 0  : in sub.spher. ==='
        Stop '=== stop in sub.spher. (l<0) ==='
! l=0
      Else If (l==0) Then
        ylm(1) = sqrt(1.0E0_dp/(4.0E0_dp*pi))
! l=1
      Else If (l==1) Then
        fac = sqrt(3.0E0_dp/(4.0E0_dp*pi))
        ylm(1) = fac*sqrt((1.0E0_dp-x*x)/2.0E0_dp)
        ylm(2) = fac*x
        ylm(3) = -ylm(1)
! l>1
      Else
        ylm(1) = 1.0E0_dp
        ylm(2) = x
        Do i = 2, l
          ylm(i+1) = ((2*i-1)*x*ylm(i)-(i-1)*ylm(i-1))/i
        End Do
        fac = 1.0E0_dp/sqrt(1.0E0_dp-x*x)
        Do m = 1, l
          lm = l + m
          ylm(lm+1) = fac*(-(l-m+1)*x*ylm(lm)+(lm-1)*ylm(l))
          If (m<l) Then
            nn = m + 1
            Do i = nn, l
              ii = l - i + nn
              ylm(ii) = fac*(-(ii-m)*x*ylm(ii)+(ii+m-2)*ylm(ii-1))
            End Do
          End If
        End Do
        fac = sqrt((2*l+1)/(4.0E0_dp*pi))
        ylm(l+1) = fac*ylm(l+1)
        Do m = 1, l
          fac = -fac/sqrt(real((l+m)*(l-m+1),kind=dp))
          lm = l + 1 + m
          ln = l + 1 - m
          qq = ylm(lm)
          ylm(lm) = fac*qq
          ylm(ln) = abs(fac)*qq
        End Do
      End If

      Return
100   Format (/, /, 3X, '==invalid argument for spher; x=', D24.16, ' ==')
    End Subroutine
