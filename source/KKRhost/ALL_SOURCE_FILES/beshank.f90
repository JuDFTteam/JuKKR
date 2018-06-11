    Subroutine beshank(hl, jl, z, lmax)
!-----------------------------------------------------------------------
!  calculates spherical bessel, hankel and neumann functions
!  for the orders lmin .le. l .le. lmax.
!  For |z| .lt. l+1 the taylor expansions of jl and nl are used.
!  For |z| .ge. l+1 the explicit expressions for hl(+), hl(-) are used.
!-----------------------------------------------------------------------
      Use mod_datatypes, Only: dp
      Implicit None
!     .. Parameters ..
      Complex (Kind=dp) :: ci
      Parameter (ci=(0.0E0_dp,1.0E0_dp))
!     ..
!     .. Scalar Arguments ..
      Complex (Kind=dp) :: z
      Integer :: lmax
!     ..
!     .. Array Arguments ..
      Complex (Kind=dp) :: hl(0:lmax), jl(0:lmax), nl(0:lmax)
!     ..
!     .. Local Scalars ..
      Complex (Kind=dp) :: termj, termn, z2, zj, zn
      Real (Kind=dp) :: rl, rn, rnm
      Integer :: l, m, n
!     ..
!     .. Intrinsic Functions ..
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

    Subroutine beshank_smallcomp(hl, jl, zval, tau, eryd, lmax)
      Use mod_datatypes, Only: dp
      Implicit None
!-----------------------------------------------------------------------
!  takes the spherical bessel etc functions stored in an array up to LMAX
!  array entries from LMAX+1 to 2*LMAX are assumed to be empty
!  these values are filled with the potential-free solution of the
!  SRA-equations
!-----------------------------------------------------------------------
      Integer :: lmax
      Complex (Kind=dp) :: hl(0:2*(lmax+1)-1), jl(0:2*(lmax+1)-1), &
        nl(0:2*(lmax+1)-1)
      Real (Kind=dp) :: cvlight
      Parameter (cvlight=274.0720442E0_dp)
      Complex (Kind=dp) :: zval
      Complex (Kind=dp) :: eryd
      Real (Kind=dp) :: tau

!       real (kind=dp) CVLIGHT
      Complex (Kind=dp) :: prefac
      Integer :: il, il2


      prefac = 1.0E0_dp/(1.0E0_dp+eryd/cvlight**2)/tau !/cvlight  !last cvlight for small component test

      il = 0
      il2 = il + lmax + 1
      nl(il2) = prefac*(zval*(-nl(il+1)))
      jl(il2) = prefac*(zval*(-jl(il+1)))
!       HL(IL2)=JL(IL2)+ CI*NL(IL2)
      hl(il2) = prefac*(zval*(-hl(il+1)))
!       write(*,'(5000E)') tau,HL(IL2),JL(IL2)+ (0.0D0,1.0D0)*NL(IL2)
!       write(*,'(5000E)') tau,HL(0),JL(0)+ (0.0D0,1.0D0)*NL(0)

      prefac = 1.0E0_dp/(1.0E0_dp+eryd/cvlight**2)/tau !/cvlight !last cvlight for small component test

      Do il = 1, lmax
        il2 = il + lmax + 1
        nl(il2) = prefac*(zval*nl(il-1)-(il+1)*nl(il))
        jl(il2) = prefac*(zval*jl(il-1)-(il+1)*jl(il))
!         HL(IL2)=JL(IL2)+ CI*NL(IL2)
        hl(il2) = prefac*(zval*hl(il-1)-(il+1)*hl(il))
!         HL(IL2)=PREFAC * ( ZVAL * HL(IL-1)-(IL+1)*HL(IL) )
!         write(*,'(5000E)') tau,HL(IL2),JL(IL2)+ (0.0D0,1.0D0)*NL(IL2)
      End Do

    End Subroutine
