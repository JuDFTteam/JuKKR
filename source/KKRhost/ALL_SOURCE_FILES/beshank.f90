SUBROUTINE beshank(hl,jl,z,lmax)
!-----------------------------------------------------------------------
!  calculates spherical bessel, hankel and neumann functions
!  for the orders lmin .le. l .le. lmax.
!  For |z| .lt. l+1 the taylor expansions of jl and nl are used.
!  For |z| .ge. l+1 the explicit expressions for hl(+), hl(-) are used.
!-----------------------------------------------------------------------
implicit none
!     .. Parameters ..
DOUBLE COMPLEX ci
PARAMETER (ci= (0.0D0,1.0D0))
!     ..
!     .. Scalar Arguments ..
DOUBLE COMPLEX z
INTEGER :: lmax
!     ..
!     .. Array Arguments ..
DOUBLE COMPLEX hl(0:lmax),jl(0:lmax),nl(0:lmax)
!     ..
!     .. Local Scalars ..
DOUBLE COMPLEX termj,termn,z2,zj,zn
DOUBLE PRECISION :: rl,rn,rnm
INTEGER :: l,m,n
!     ..
!     .. Intrinsic Functions ..
INTRINSIC CDABS,EXP
!     ..
zj = 1.d0
zn = 1.d0
z2 = z*z
IF (CDABS(z) < lmax+1.d0) THEN
  DO  l = 0,lmax
    rl = l + l
    termj = -0.5D0/ (rl+3.d0)*z2
    termn = 0.5D0/ (rl-1.d0)*z2
    jl(l) = 1.d0
    nl(l) = 1.d0
    DO  n = 2,25
      jl(l) = jl(l) + termj
      nl(l) = nl(l) + termn
      rn = n + n
      termj = -termj/ (rl+rn+1.d0)/rn*z2
      termn = termn/ (rl-rn+1.d0)/rn*z2
    END DO
    jl(l) = jl(l)*zj
    nl(l) = -nl(l)*zn/z
    hl(l) = jl(l) + nl(l)*ci
    
    zj = zj*z/ (rl+3.d0)
    zn = zn/z* (rl+1.d0)
  END DO
END IF

DO  l = 0,lmax
  IF (CDABS(z) >= l+1.d0) THEN
    hl(l) = 0.d0
    nl(l) = 0.d0
    rnm = 1.d0
    DO  m = 0,l
      hl(l) = hl(l) + rnm/ (-ci* (z+z))**m
      nl(l) = nl(l) + rnm/ (ci* (z+z))**m
      rnm = rnm* (l*l+l-m*m-m)/ (m+1.d0)
    END DO
    hl(l) = hl(l)* (-ci)**l*EXP(ci*z)/ (ci*z)
    nl(l) = nl(l)*ci**l*EXP(-ci*z)/ (-ci*z)
    jl(l) = (hl(l)+nl(l))*0.5D0
    nl(l) = (hl(l)-jl(l))/ci
  END IF
END DO

RETURN

END SUBROUTINE

SUBROUTINE beshank_smallcomp(hl,jl,zval,tau,eryd,lmax)
IMPLICIT NONE
!-----------------------------------------------------------------------
!  takes the spherical bessel etc functions stored in an array up to LMAX
!  array entries from LMAX+1 to 2*LMAX are assumed to be empty
!  these values are filled with the potential-free solution of the
!  SRA-equations
!-----------------------------------------------------------------------
DOUBLE COMPLEX hl(0:2*(lmax+1)-1), jl(0:2*(lmax+1)-1),  &
    nl(0:2*(lmax+1)-1)
DOUBLE PRECISION :: cvlight
PARAMETER (cvlight=274.0720442D0)
DOUBLE COMPLEX zval
DOUBLE COMPLEX eryd
DOUBLE PRECISION :: tau
INTEGER :: lmax

!       DOUBLE PRECISION CVLIGHT
DOUBLE COMPLEX prefac
INTEGER :: il,il2


prefac = 1.0D0 / (1.0D0+eryd/cvlight**2) / tau !/cvlight  !last cvlight for small component test

il=0
il2=il+lmax+1
nl(il2)=prefac * (zval* (-nl(il+1)) )
jl(il2)=prefac * (zval* (-jl(il+1)) )
!       HL(IL2)=JL(IL2)+ CI*NL(IL2)
hl(il2)=prefac * (zval* (-hl(il+1)) )
!       write(*,'(5000E)') tau,HL(IL2),JL(IL2)+ (0.0D0,1.0D0)*NL(IL2)
!       write(*,'(5000E)') tau,HL(0),JL(0)+ (0.0D0,1.0D0)*NL(0)

prefac = 1.0D0 / (1.0D0+eryd/cvlight**2) / tau !/cvlight !last cvlight for small component test

DO il=1,lmax
  il2=il+lmax+1
  nl(il2)=prefac * ( zval * nl(il-1)-(il+1)*nl(il) )
  jl(il2)=prefac * ( zval * jl(il-1)-(il+1)*jl(il) )
!         HL(IL2)=JL(IL2)+ CI*NL(IL2)
  hl(il2)=prefac * ( zval * hl(il-1)-(il+1)*hl(il) )
!         HL(IL2)=PREFAC * ( ZVAL * HL(IL-1)-(IL+1)*HL(IL) )
!         write(*,'(5000E)') tau,HL(IL2),JL(IL2)+ (0.0D0,1.0D0)*NL(IL2)
END DO

END SUBROUTINE beshank_smallcomp

