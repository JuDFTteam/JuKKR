SUBROUTINE beshan(hl,jl,nl,z,lmax)
!-----------------------------------------------------------------------
!  calculates spherical bessel, hankel and neumann functions
!  for the orders l .le. lmax.
!  For |z| .lt. 1 the taylor expansions of jl and nl are used.
!  For |z| .ge. 1 the explicit expressions for hl(+), hl(-) are used.

!                            R. Zeller   Jan. 1990
!-----------------------------------------------------------------------
!     .. Parameters ..
      DOUBLE COMPLEX CI
      PARAMETER (CI= (0.0D0,1.0D0))
!..
!.. Scalar Arguments ..
      DOUBLE COMPLEX Z
      INTEGER LMAX
!..
!.. Array Arguments ..
      DOUBLE COMPLEX HL(0:LMAX),JL(0:LMAX),NL(0:LMAX)
!..
!.. Local Scalars ..
      DOUBLE COMPLEX TERMJ,TERMN,Z2,ZJ,ZN
      DOUBLE PRECISION RL,RN,RNM
      INTEGER L,M,N
!..
!.. Intrinsic Functions ..
      INTRINSIC ABS,EXP

!     ..
zj = 1.d0
zn = 1.d0
z2 = z*z
IF (ABS(z) < lmax+1.d0) THEN
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
  IF (ABS(z) >= l+1.d0) THEN
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

END SUBROUTINE beshan
