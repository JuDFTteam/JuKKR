function cdnlzdz(l, z, mode)
!   ********************************************************************
!   *                                                                  *
!   *     d n(L,Z) / dz    analytically                                *
!   *                                                                  *
!   ********************************************************************
  implicit none

! Dummy arguments
  integer :: l, mode
  complex *16 :: z
  complex *16 :: cdnlzdz

! Local variables
  complex *16 :: cnlz

  if (mode==1) then

    if (l==0) then

      cdnlzdz = l*cnlz(l, z)/z - cnlz(l+1, z)
    else
      cdnlzdz = (l*cnlz(l-1,z)-(l+1)*cnlz(l+1,z))/dble(2*l+1)
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
end function
