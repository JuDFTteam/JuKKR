! **********************************************************************
    Subroutine ymy(v1, v2, v3, r, ylm, lmax)
      Use mod_datatypes, Only: dp
! **********************************************************************
!    this subroutine calculates real spherical harmonics with the
!     normalization : <y|y> =1
!    returns also r = length of vector v

!     generate the complex spherical harmonics for the vector v
!     using a stable upward recursion in l.  (see notes
!     by m. weinert.)
!                                  m.weinert  1982

!     converted to real spherical harmonics .
!                                  b.drittler 1987
!-----------------------------------------------------------------------

      Implicit None
!.. Parameters ..
      Real (Kind=dp) :: szero
      Parameter (szero=1.0E-20_dp)
!..
!.. Scalar Arguments ..
      Real (Kind=dp), Intent (Out) :: r
      Real (Kind=dp), Intent (In) :: v1, v2, v3
      Integer, Intent (In) :: lmax
!..
!.. Array Arguments ..
      Real (Kind=dp), Intent (Out) :: ylm((2*lmax+1)**2)
!..
!.. Local Scalars ..
      Real (Kind=dp) :: a, cd, cph, cth, fac, fpi, pi, rtwo, sgm, sph, sth, t, &
        xy, xyz
      Integer :: i, l, m
!..
!.. Local Arrays ..
      Real (Kind=dp) :: c(0:lmax), p(0:lmax, 0:lmax), s(0:lmax)
!..
!.. Intrinsic Functions ..
      Intrinsic :: atan, sqrt
!..
!.. External Subroutines ..
      External :: rcstop
!     ..
      pi = 4.E0_dp*atan(1.E0_dp)
      fpi = 4.E0_dp*pi
      rtwo = sqrt(2.0E0_dp)

!--->    calculate sin and cos of theta and phi

      xy = v1**2 + v2**2
      xyz = xy + v3**2

      r = sqrt(xyz)
      If (xyz<=0.0E0_dp) Then
        Write (*, *) xyz
        Call rcstop('ylm=0   ')

      Else

        If (xy>szero*xyz) Then
          xy = sqrt(xy)
          xyz = sqrt(xyz)
          cth = v3/xyz
          sth = xy/xyz
          cph = v1/xy
          sph = v2/xy

        Else

          sth = 0.0E0_dp
          cth = 1.0E0_dp
          If (v3<0) cth = -1.0E0_dp
          cph = 1.0E0_dp
          sph = 0.0E0_dp
        End If

!--->    generate associated legendre functions for m.ge.0
!        loop over m values

        fac = 1.0E0_dp
        Do m = 0, lmax - 1
          fac = -(2*m-1)*fac
          p(m, m) = fac
          p(m+1, m) = (2*m+1)*cth*fac

!--->    recurse upward in l

          Do l = m + 2, lmax
            p(l, m) = ((2*l-1)*cth*p(l-1,m)-(l+m-1)*p(l-2,m))/(l-m)
          End Do
          fac = fac*sth
        End Do
        p(lmax, lmax) = -(2*lmax-1)*fac

!--->    determine sin and cos of phi

        s(0) = 0.0E0_dp
        s(1) = sph
        c(0) = 1.0E0_dp
        c(1) = cph
        Do m = 2, lmax
          s(m) = 2*cph*s(m-1) - s(m-2)
          c(m) = 2*cph*c(m-1) - c(m-2)
        End Do

!--->    multiply in the normalization factors

        i = 0
        Do l = 0, lmax
          i = i + l + 1
          a = sqrt((2*l+1)/fpi)
          cd = 1
          ylm(i) = a*p(l, 0)
          sgm = -rtwo
          Do m = 1, l
            t = (l+1-m)*(l+m)
            cd = cd/t
            t = a*sqrt(cd)
            ylm(i+m) = sgm*t*p(l, m)*c(m)
            ylm(i-m) = sgm*t*p(l, m)*s(m)
            sgm = -sgm
          End Do
          i = i + l
        End Do

      End If

      Return

    End Subroutine
