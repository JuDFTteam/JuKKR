module mod_ioben
  
  private
  public :: ioben

contains

  !-------------------------------------------------------------------------------
  !> Summary: Return next upper or equal integer
  !> Author: Rudolf Berrendorf
  !> Date: July 1992
  !> Category: KKRhost, numerical-tools
  !> Deprecated: True ! This needs to be set to True for deprecated subroutines
  !>
  !>                         --   --
  !> Calculates the function |  r  |  (next upper or equal integer)
  !>                         |     |
  !>
  !> Descrition of input parameters:
  !>
  !> r : real number to look for
  !>
  !> last update: February 1994
  !> 
  !> @note seems to be unused @endnote
  !-------------------------------------------------------------------------------
  integer function ioben(r)

    use :: mod_datatypes, only: dp
    implicit none

    ! .. Scalar Arguments ..
    real (kind=dp) :: r
    ! ..
    ! .. Intrinsic Functions ..
    intrinsic :: abs, int, nint
    ! ..

    ! .. Parameters ..
    real (kind=dp), parameter :: eps=1e-6_dp
    ! ..

    if ((nint(r)-r)<eps) then
      if (abs(nint(r)-r)<eps) then
        ioben = nint(r)
      else
        ioben = nint(r+1.0_dp)
      end if
    else
      if (abs(int(r)-r)<eps) then
        ioben = int(r)
      else
        ioben = int(r+1.0_dp)
      end if
    end if

    return
  end function ioben

end module mod_ioben
