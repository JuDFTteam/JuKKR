subroutine beshank(hl, jl, z, lmax)
!-----------------------------------------------------------------------
!  calculates spherical bessel, hankel and neumann functions
!  for the orders lmin .le. l .le. lmax.
!  For |z| .lt. l+1 the taylor expansions of jl and nl are used.
!  For |z| .ge. l+1 the explicit expressions for hl(+), hl(-) are used.
!-----------------------------------------------------------------------
  implicit none
!     .. Parameters ..
  double complex :: ci
  parameter (ci=(0.0d0,1.0d0))
!     ..
!     .. Scalar Arguments ..
  double complex :: z
  integer :: lmax
!     ..
!     .. Array Arguments ..
  double complex :: hl(0:lmax), jl(0:lmax), nl(0:lmax)
!     ..
!     .. Local Scalars ..
  double complex :: termj, termn, z2, zj, zn
  double precision :: rl, rn, rnm
  integer :: l, m, n
!     ..
!     .. Intrinsic Functions ..
  intrinsic :: cdabs, exp
!     ..
  zj = 1.d0
  zn = 1.d0
  z2 = z*z
  if (cdabs(z)<lmax+1.d0) then
    do l = 0, lmax
      rl = l + l
      termj = -0.5d0/(rl+3.d0)*z2
      termn = 0.5d0/(rl-1.d0)*z2
      jl(l) = 1.d0
      nl(l) = 1.d0
      do n = 2, 25
        jl(l) = jl(l) + termj
        nl(l) = nl(l) + termn
        rn = n + n
        termj = -termj/(rl+rn+1.d0)/rn*z2
        termn = termn/(rl-rn+1.d0)/rn*z2
      end do
      jl(l) = jl(l)*zj
      nl(l) = -nl(l)*zn/z
      hl(l) = jl(l) + nl(l)*ci

      zj = zj*z/(rl+3.d0)
      zn = zn/z*(rl+1.d0)
    end do
  end if

  do l = 0, lmax
    if (cdabs(z)>=l+1.d0) then
      hl(l) = 0.d0
      nl(l) = 0.d0
      rnm = 1.d0
      do m = 0, l
        hl(l) = hl(l) + rnm/(-ci*(z+z))**m
        nl(l) = nl(l) + rnm/(ci*(z+z))**m
        rnm = rnm*(l*l+l-m*m-m)/(m+1.d0)
      end do
      hl(l) = hl(l)*(-ci)**l*exp(ci*z)/(ci*z)
      nl(l) = nl(l)*ci**l*exp(-ci*z)/(-ci*z)
      jl(l) = (hl(l)+nl(l))*0.5d0
      nl(l) = (hl(l)-jl(l))/ci
    end if
  end do

  return

end subroutine

subroutine beshank_smallcomp(hl, jl, zval, tau, eryd, lmax)
  implicit none
!-----------------------------------------------------------------------
!  takes the spherical bessel etc functions stored in an array up to LMAX
!  array entries from LMAX+1 to 2*LMAX are assumed to be empty
!  these values are filled with the potential-free solution of the
!  SRA-equations
!-----------------------------------------------------------------------
  double complex :: hl(0:2*(lmax+1)-1), jl(0:2*(lmax+1)-1), nl(0:2*(lmax+1)-1)
  double precision :: cvlight
  parameter (cvlight=274.0720442d0)
  double complex :: zval
  double complex :: eryd
  double precision :: tau
  integer :: lmax

!       DOUBLE PRECISION CVLIGHT
  double complex :: prefac
  integer :: il, il2


  prefac = 1.0d0/(1.0d0+eryd/cvlight**2)/tau !/cvlight  !last cvlight for small component test

  il = 0
  il2 = il + lmax + 1
  nl(il2) = prefac*(zval*(-nl(il+1)))
  jl(il2) = prefac*(zval*(-jl(il+1)))
!       HL(IL2)=JL(IL2)+ CI*NL(IL2)
  hl(il2) = prefac*(zval*(-hl(il+1)))
!       write(*,'(5000E)') tau,HL(IL2),JL(IL2)+ (0.0D0,1.0D0)*NL(IL2)
!       write(*,'(5000E)') tau,HL(0),JL(0)+ (0.0D0,1.0D0)*NL(0)

  prefac = 1.0d0/(1.0d0+eryd/cvlight**2)/tau !/cvlight !last cvlight for small component test

  do il = 1, lmax
    il2 = il + lmax + 1
    nl(il2) = prefac*(zval*nl(il-1)-(il+1)*nl(il))
    jl(il2) = prefac*(zval*jl(il-1)-(il+1)*jl(il))
!         HL(IL2)=JL(IL2)+ CI*NL(IL2)
    hl(il2) = prefac*(zval*hl(il-1)-(il+1)*hl(il))
!         HL(IL2)=PREFAC * ( ZVAL * HL(IL-1)-(IL+1)*HL(IL) )
!         write(*,'(5000E)') tau,HL(IL2),JL(IL2)+ (0.0D0,1.0D0)*NL(IL2)
  end do

end subroutine

