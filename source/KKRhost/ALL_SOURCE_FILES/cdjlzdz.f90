function cdjlzdz(l, z, mode)
  ! ********************************************************************
  ! *                                                                  *
  ! *     d j(L,Z) / dz    analytically                                *
  ! *                                                                  *
  ! ********************************************************************
  use :: mod_datatypes, only: dp
  implicit none

  ! Dummy arguments
  integer :: l, mode
  complex (kind=dp) :: z
  complex (kind=dp) :: cdjlzdz

  ! Local variables
  complex (kind=dp) :: cjlz

  if (mode==1) then

    if (l==0) then

      cdjlzdz = l*cjlz(l, z)/z - cjlz(l+1, z)
    else
      cdjlzdz = (l*cjlz(l-1,z)-(l+1)*cjlz(l+1,z))/real(2*l+1, kind=dp)
      return
    end if
  else if (mode==2) then

    if (l==0) then
      cdjlzdz = l*cjlz(l, z)/z - cjlz(l+1, z)
    else
      cdjlzdz = cjlz(l-1, z) - (l+1)*cjlz(l, z)/z
      return
    end if
  else
    cdjlzdz = l*cjlz(l, z)/z - cjlz(l+1, z)
  end if
end function cdjlzdz
