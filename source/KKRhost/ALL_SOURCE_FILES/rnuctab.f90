!------------------------------------------------------------------------------------
!> Summary: Stops the program if something goes wrong in the dirac solver
!> Author: 
!> Stops the program if something goes wrong in the dirac solver
!------------------------------------------------------------------------------------
module mod_rnuctab
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Stops the program if something goes wrong in the dirac solver
  !> Author: 
  !> Category: dirac, sanity-check, KKRhost
  !> Deprecated: False 
  !> Stops the program if something goes wrong in the dirac solver
  !-------------------------------------------------------------------------------
  function rnuctab(z)
    real (kind=dp) :: rnuctab
    integer, intent (inout) :: z

    z = 0
    rnuctab = 0e0_dp
    stop ' < RNUCTAB > : NUCLEUS <> 0 not implemented '
  end function rnuctab

end module mod_rnuctab
