module mod_dclock

contains

! **********************************************************************
function dclock()
  use :: mod_datatypes, only: dp
  real (kind=dp) :: dclock
  ! **********************************************************************
  ! .. External Functions ..
  real (kind=dp) :: etime, tarry(2)
  external :: etime

  ! ..
  dclock = real(etime(tarry), kind=dp)
  return
end function dclock

end module mod_dclock
