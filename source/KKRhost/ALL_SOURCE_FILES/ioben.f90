module mod_ioben
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: 
  !> Author: 
  !> Category: KKRhost, 
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> 
  !-------------------------------------------------------------------------------
  integer function ioben(r)
    ! -----------------------------------------------------------------------

    ! --   --
    ! Calculates the function |  r  |  (next upper or equal integer)
    ! |     |

    ! Descrition of input parameters:

    ! r : real number to look for

    ! Rudolf Berrendorf, July 1992
    ! last update: February 1994
    ! -----------------------------------------------------------------------

    implicit none

    ! .. Scalar Arguments ..

    real (kind=dp) :: r
    ! ..
    ! .. Intrinsic Functions ..
    intrinsic :: abs, int, nint
    ! ..

    ! .. Parameters ..
    real (kind=dp) :: eps
    parameter (eps=1e-6_dp)
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
