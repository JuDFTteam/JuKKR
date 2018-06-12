subroutine beshan(hl, jl, nl, z, lmax)
  ! -----------------------------------------------------------------------
  ! calculates spherical bessel, hankel and neumann functions
  ! for the orders l .le. lmax.
  ! For |z| .lt. 1 the taylor expansions of jl and nl are used.
  ! For |z| .ge. 1 the explicit expressions for hl(+), hl(-) are used.

  ! R. Zeller   Jan. 1990
  ! -----------------------------------------------------------------------
  use :: mod_datatypes, only: dp
  implicit none
  ! .. Parameters ..
  complex (kind=dp) :: ci
  parameter (ci=(0.0e0_dp,1.0e0_dp))
  ! ..
  ! .. Scalar Arguments ..
  complex (kind=dp) :: z
  integer :: lmax
  ! ..
  ! .. Array Arguments ..
  complex (kind=dp) :: hl(0:lmax), jl(0:lmax), nl(0:lmax)
  ! ..
  ! .. Local Scalars ..
  complex (kind=dp) :: termj, termn, z2, zj, zn
  real (kind=dp) :: rl, rn, rnm
  integer :: l, m, n
  ! ..
  ! .. Intrinsic Functions ..
  intrinsic :: abs, exp

  ! ..
  zj = 1.e0_dp
  zn = 1.e0_dp
  z2 = z*z
  if (abs(z)<lmax+1.e0_dp) then
    do l = 0, lmax
      rl = l + l
      termj = -0.5e0_dp/(rl+3.e0_dp)*z2
      termn = 0.5e0_dp/(rl-1.e0_dp)*z2
      jl(l) = 1.e0_dp
      nl(l) = 1.e0_dp
      do n = 2, 25
        jl(l) = jl(l) + termj
        nl(l) = nl(l) + termn
        rn = n + n
        termj = -termj/(rl+rn+1.e0_dp)/rn*z2
        termn = termn/(rl-rn+1.e0_dp)/rn*z2
      end do
      jl(l) = jl(l)*zj
      nl(l) = -nl(l)*zn/z
      hl(l) = jl(l) + nl(l)*ci

      zj = zj*z/(rl+3.e0_dp)
      zn = zn/z*(rl+1.e0_dp)
    end do
  end if

  do l = 0, lmax
    if (abs(z)>=l+1.e0_dp) then
      hl(l) = 0.e0_dp
      nl(l) = 0.e0_dp
      rnm = 1.e0_dp
      do m = 0, l
        hl(l) = hl(l) + rnm/(-ci*(z+z))**m
        nl(l) = nl(l) + rnm/(ci*(z+z))**m
        rnm = rnm*(l*l+l-m*m-m)/(m+1.e0_dp)
      end do
      hl(l) = hl(l)*(-ci)**l*exp(ci*z)/(ci*z)
      nl(l) = nl(l)*ci**l*exp(-ci*z)/(-ci*z)
      jl(l) = (hl(l)+nl(l))*0.5e0_dp
      nl(l) = (hl(l)-jl(l))/ci
    end if
  end do

  return

end subroutine beshan
