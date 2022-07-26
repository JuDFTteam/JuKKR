!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------
!> Summary: Spherical Bessel and Hanke functions for large and small components
!>
!> @note
!> Can probably be merged with beshan module
!> @endnote
!-------------------------------------------------------------------------------
module mod_beshank

  private
  public :: beshank, beshank_smallcomp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Spherical Bessel and Hanke functions
  !> Author: 
  !> Category: KKRhost, special-functions, single-site
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Calculates spherical Bessel, Hankel and Neumann (not returned) functions
  !> for the orders lmin .le. l .le. lmax.
  !> For |z| < l+1 the Taylor expansions of jl and nl are used.
  !> For |z| >= l+1 the explicit expressions for hl(+), hl(-) are used.
  !-------------------------------------------------------------------------------
  subroutine beshank(hl, jl, z, lmax)
    use :: mod_datatypes, only: dp
    use :: mod_constants, only: ci
    implicit none

    ! inputs
    complex (kind=dp), intent (in) :: z !! complex energy E+i*eta = Re(z)+i*Im(z)
    integer, intent (in) :: lmax   !! energy cutoff

    ! outputs
    complex (kind=dp), intent (out) :: hl(0:lmax), jl(0:lmax) !! spherical Hankel and Bessel up to lmax

    ! locals
    complex (kind=dp) :: nl(0:lmax) !! Neumannc function
    complex (kind=dp) :: termj, termn, z2, zj, zn
    real (kind=dp) :: rl, rn, rnm
    integer :: l, m, n


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
          ! rn = real(n + n, kind=dp)
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

  end subroutine beshank


  !-------------------------------------------------------------------------------
  !> Summary: Spherical Bessel and Hanke functions for small component
  !> Author: 
  !> Category: KKRhost, special-functions, single-site
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Takes the spherical bessel etc functions stored in an array up to LMAX
  !> array entries from LMAX+1 to 2*LMAX are assumed to be empty
  !> these values are filled with the potential-free solution of the
  !> SRA-equations
  !-------------------------------------------------------------------------------
  subroutine beshank_smallcomp(hl, jl, zval, tau, eryd, lmax)
    use :: mod_datatypes, only: dp
    implicit none
    integer, intent (in) :: lmax
    complex (kind=dp), intent (in) :: zval
    complex (kind=dp), intent (in) :: eryd
    real (kind=dp), intent (in) :: tau
    complex (kind=dp), intent (inout) :: hl(0:2*(lmax+1)-1)
    complex (kind=dp), intent (inout) :: jl(0:2*(lmax+1)-1)
    ! complex (kind=dp), intent(inout) :: nl(0:2*(lmax+1)-1)

    ! local
    real (kind=dp), parameter :: cvlight = 274.0720442e0_dp
    complex (kind=dp) :: prefac
    integer :: il, il2


    prefac = 1.0e0_dp/(1.0e0_dp+eryd/cvlight**2)/tau ! last cvlight for small component test

    il = 0
    il2 = il + lmax + 1
    ! nl(il2) = prefac*(zval*(-nl(il+1)))
    jl(il2) = prefac*(zval*(-jl(il+1)))
    ! HL(IL2)=JL(IL2)+ CI*NL(IL2)
    hl(il2) = prefac*(zval*(-hl(il+1)))

    prefac = 1.0e0_dp/(1.0e0_dp+eryd/cvlight**2)/tau ! last cvlight for small component test

    do il = 1, lmax
      il2 = il + lmax + 1
      ! nl(il2) = prefac*(zval*nl(il-1)-(il+1)*nl(il))
      jl(il2) = prefac*(zval*jl(il-1)-(il+1)*jl(il))
      ! HL(IL2)=JL(IL2)+ CI*NL(IL2)
      hl(il2) = prefac*(zval*hl(il-1)-(il+1)*hl(il))
    end do

  end subroutine beshank_smallcomp

end module mod_beshank
