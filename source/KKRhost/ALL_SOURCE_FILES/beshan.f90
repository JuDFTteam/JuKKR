    Subroutine beshan(hl, jl, nl, z, lmax)
      Use mod_datatypes, Only: dp
!-----------------------------------------------------------------------
!  calculates spherical bessel, hankel and neumann functions
!  for the orders l .le. lmax.
!  For |z| .lt. 1 the taylor expansions of jl and nl are used.
!  For |z| .ge. 1 the explicit expressions for hl(+), hl(-) are used.

!                            R. Zeller   Jan. 1990
!-----------------------------------------------------------------------
!     .. Parameters ..
      Complex (Kind=dp) :: ci
      Parameter (ci=(0.0E0_dp,1.0E0_dp))
!..
!.. Scalar Arguments ..
      Complex (Kind=dp) :: z
      Integer :: lmax
!..
!.. Array Arguments ..
      Complex (Kind=dp) :: hl(0:lmax), jl(0:lmax), nl(0:lmax)
!..
!.. Local Scalars ..
      Complex (Kind=dp) :: termj, termn, z2, zj, zn
      Real (Kind=dp) :: rl, rn, rnm
      Integer :: l, m, n
!..
!.. Intrinsic Functions ..
      Intrinsic :: abs, exp

!     ..
      zj = 1.E0_dp
      zn = 1.E0_dp
      z2 = z*z
      If (abs(z)<lmax+1.E0_dp) Then
        Do l = 0, lmax
          rl = l + l
          termj = -0.5E0_dp/(rl+3.E0_dp)*z2
          termn = 0.5E0_dp/(rl-1.E0_dp)*z2
          jl(l) = 1.E0_dp
          nl(l) = 1.E0_dp
          Do n = 2, 25
            jl(l) = jl(l) + termj
            nl(l) = nl(l) + termn
            rn = n + n
            termj = -termj/(rl+rn+1.E0_dp)/rn*z2
            termn = termn/(rl-rn+1.E0_dp)/rn*z2
          End Do
          jl(l) = jl(l)*zj
          nl(l) = -nl(l)*zn/z
          hl(l) = jl(l) + nl(l)*ci

          zj = zj*z/(rl+3.E0_dp)
          zn = zn/z*(rl+1.E0_dp)
        End Do
      End If

      Do l = 0, lmax
        If (abs(z)>=l+1.E0_dp) Then
          hl(l) = 0.E0_dp
          nl(l) = 0.E0_dp
          rnm = 1.E0_dp
          Do m = 0, l
            hl(l) = hl(l) + rnm/(-ci*(z+z))**m
            nl(l) = nl(l) + rnm/(ci*(z+z))**m
            rnm = rnm*(l*l+l-m*m-m)/(m+1.E0_dp)
          End Do
          hl(l) = hl(l)*(-ci)**l*exp(ci*z)/(ci*z)
          nl(l) = nl(l)*ci**l*exp(-ci*z)/(-ci*z)
          jl(l) = (hl(l)+nl(l))*0.5E0_dp
          nl(l) = (hl(l)-jl(l))/ci
        End If
      End Do

      Return

    End Subroutine
