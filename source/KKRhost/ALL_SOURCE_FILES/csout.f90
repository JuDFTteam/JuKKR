module mod_csout
  
  private
  public :: csout

contains

  !-------------------------------------------------------------------------------
  !> Summary: Outward integration of multiple functions with ext. 3-point Simpson
  !> Author: B. Drittler, M. Ogura
  !> Date: 1989
  !> Category: KKRhost, radial-grid, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> This subroutine does an outwards integration of llmax
  !> functions f with an extended 3-point-simpson:
  !>
  !> ir
  !> fint(ll,i) = { f(ll,i') di'
  !> irmin
  !>
  !> The starting value for this integration at irmin+1 is determined by
  !> a 4 point lagrangian integration  , coefficients given by
  !> m. abramowitz and i.a. stegun, handbook of mathematical functions,
  !> nbs applied mathematics series 55 (1968)
  !>
  !> @warning In case of radial integration :
  !> the weights drdi have to be multiplied before calling this
  !> subroutine.
  !> @endwarning
  !>
  !> B. Drittler Mar. 1989
  !>
  !> Modified for functions with kinks - at each kink the integration
  !> is restarted
  !>
  !> @warning It is supposed that irmin + 3 is less than imt! @endwarning
  !>
  !> B. Drittler july 1989
  !> Modified by M. Ogura, June 2015
  !>
  !> @note
  !> Maybe this is a duplication of `csimpk`?
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine csout(f, fint, lmmsqd, irmind, irmd, irmin, ipan, ircut)

    use :: mod_datatypes, only: dp
    implicit none

    ! .. Parameters ..
    real (kind=dp) :: a1, a2, a3
    parameter (a1=5.e0_dp/12.e0_dp, a2=8.e0_dp/12.e0_dp, a3=-1.e0_dp/12.e0_dp)
    ! ..
    ! .. Scalar Arguments ..
    integer :: ipan, irmd, irmind, lmmsqd, irmin
    ! ..
    ! .. Array Arguments ..
    complex (kind=dp) :: f(lmmsqd, irmind:irmd), fint(lmmsqd, irmind:irmd)
    integer :: ircut(0:ipan)
    ! ..
    ! .. Local Scalars ..
    integer :: i, ien, ip, ist, ll

    ! ---> loop over kinks

    do ip = 1, ipan
      ien = ircut(ip)
      ist = ircut(ip-1) + 1

      if (ip==1) then
        ist = irmin
        do ll = 1, lmmsqd
          fint(ll, ist) = 0.0e0_dp
        end do

      else
        do ll = 1, lmmsqd
          fint(ll, ist) = fint(ll, ist-1)
        end do
      end if


      ! ---> calculate fint with an extended 3-point-simpson

      do i = ist, ien - 2, 2
        do ll = 1, lmmsqd
          fint(ll, i+1) = fint(ll, i) + f(ll, i)*a1 + f(ll, i+1)*a2 + f(ll, i+2)*a3
          fint(ll, i+2) = fint(ll, i+1) + f(ll, i)*a3 + f(ll, i+1)*a2 + f(ll, i+2)*a1
        end do
      end do
      if (mod(ien-ist,2)==1) then
        do ll = 1, lmmsqd
          fint(ll, ien) = fint(ll, ien-1) + f(ll, ien)*a1 + f(ll, ien-1)*a2 + f(ll, ien-2)*a3
        end do
      end if
    end do

  end subroutine csout

end module mod_csout
