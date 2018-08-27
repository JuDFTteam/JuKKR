module mod_rnuctab

contains

real (kind=dp) function rnuctab(z)
  use :: mod_datatypes, only: dp
  integer, intent (inout) :: z

  z = 0
  rnuctab = 0e0_dp
  stop ' < RNUCTAB > : NUCLEUS <> 0 not implemented '
end function rnuctab

end module mod_rnuctab
