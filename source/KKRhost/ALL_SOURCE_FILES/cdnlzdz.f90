function cdnlzdz(l, z, mode)
  ! ********************************************************************
  ! *                                                                  *
  ! *     d n(L,Z) / dz    analytically                                *
  ! *                                                                  *
  ! ********************************************************************
  use :: mod_datatypes, only: dp
  implicit none

  ! Dummy arguments
  integer :: l, mode
  complex (kind=dp) :: z
  complex (kind=dp) :: cdnlzdz

  ! Local variables
  complex (kind=dp) :: cnlz

  if (mode==1) then

    if (l==0) then

      cdnlzdz = l*cnlz(l, z)/z - cnlz(l+1, z)
    else
      cdnlzdz = (l*cnlz(l-1,z)-(l+1)*cnlz(l+1,z))/real(2*l+1, kind=dp)
      return
    end if
  else if (mode==2) then

    if (l==0) then
      cdnlzdz = l*cnlz(l, z)/z - cnlz(l+1, z)
    else
      cdnlzdz = cnlz(l-1, z) - (l+1)*cnlz(l, z)/z
      return
    end if
  else
    cdnlzdz = l*cnlz(l, z)/z - cnlz(l+1, z)
  end if
end function cdnlzdz
