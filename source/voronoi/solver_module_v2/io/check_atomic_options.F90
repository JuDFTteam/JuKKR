  subroutine check_atomic_options(itc)
! Generates text describing the atomic options set by inpsusc.dat
! Checks consistency
  use global

  implicit none

  integer(kind=i4b), intent(in) :: itc
! --------------------------------------
  integer(kind=i4b) :: ia, il, ib

  nasusc2 = 0
! **************
  do ia=1,nasusc
! **************
    if (itc == 1) write(*,'("************************************************************")')
    if (itc == 1) write(*,'("  ia=",i4,"  iasusc=",i4)') ia, iasusc(ia)
!   SOC option for current atom
    if (lsoc) then 
      if (isoc(ia) == 1) then
        if (itc == 1) write(*,'("  SOC full: L.S")')
        if (itc == 1) write(*,'("  SOC applied with scaling=",es12.4)') socscaling(ia)
      else if (isoc(ia) == 2) then
        if (itc == 1) write(*,'("  SOC longitudinal: (L.u).(S.u)")')
        if (itc == 1) write(*,'("  SOC applied with scaling=",es12.4)') socscaling(ia)
      else if (isoc(ia) == 3) then
        if (itc == 1) write(*,'("  SOC transverse: L.S - (L.u).(S.u)")')
        if (itc == 1) write(*,'("  SOC applied with scaling=",es12.4)') socscaling(ia)
      else if (isoc(ia) == 0) then
        if (itc == 1) write(*,'("  No SOC applied")')
      else
        stop 'Unknown SOC option'
      end if
    end if
!   External B field for current atom
    if (lbfield) then 
      if (ibfield(ia) == 1) then
        if (itc == 1) write(*,'("  B field in B.S form")')
        if (itc == 1) write(*,'("  B field applied=",3es12.4)') blen(ia)*bdir(:,ia)
      else if (ibfield(ia) == 2) then
        if (itc == 1) write(*,'("  B field in B.L form")')
        if (itc == 1) write(*,'("  B field applied=",3es12.4)') blen(ia)*bdir(:,ia)
      else if (ibfield(ia) == 3) then
        if (itc == 1) write(*,'("  B field in B.(L+S) form")')
        if (itc == 1) write(*,'("  B field applied=",3es12.4)') blen(ia)*bdir(:,ia)
      else if (ibfield(ia) == 0) then
        if (itc == 1) write(*,'("  No B field applied")')
      else
        stop 'Unknown B field option'
      end if
    end if
!   Susceptibility options
    if (lsusc) then
      if (isusc(ia) == 1 .or. isusc(ia) == -1) then
        if (itc == 1) write(*,'("  Transverse susceptibility")')
      else if (isusc(ia) == 2 .or. isusc(ia) == -2) then
        if (itc == 1) write(*,'("  Longitudinal susceptibility")')
      else if (isusc(ia) == 3 .or. isusc(ia) == -3) then
        if (itc == 1) write(*,'("  Full susceptibility")')
      else if (isusc(ia) == 0) then
        if (itc == 1) write(*,'("  No susceptibility")')
      else
        stop 'Unknown susceptibility option'
      end if
!     Which atoms for the susceptibility
      if (isusc(ia) /= 0) then
        nasusc2 = nasusc2 + 1
        iasusc2(nasusc2) = ia
!       Hartree kernel options
        if (lkha) then
          if (ikha(ia) == 1) then
            if (itc == 1) write(*,'("  Hartree kernel")')
          else if (ikha(ia) == 0) then
            if (itc == 1) write(*,'("  No Hartree kernel")')
          else
            stop 'Unknown Hartree kernel option'
          end if
        end if
!       xc kernel options
        if (lkxc) then
          if (ikxc(ia) == 1) then
            if (itc == 1) write(*,'("  Transverse xc kernel")')
          else if (ikxc(ia) == 2) then
            if (itc == 1) write(*,'("  Longitudinal xc kernel")')
          else if (ikxc(ia) == 3) then
            if (itc == 1) write(*,'("  Full xc kernel")')
          else if (ikxc(ia) == 0) then
            if (itc == 1) write(*,'("  No xc kernel")')
          else
            stop 'Unknown xc kernel option'
          end if
        end if
      end if
    end if
! ******
  end do
! ******
  if (itc == 1) write(*,'("************************************************************")')
  if (itc == 1 .and. lsusc) write(*,'(" Number of atoms for susceptibility calculation:",i4)') nasusc2
  if (itc == 1 .and. lsusc) write(*,'(" iasusc2=",1000i4)') iasusc2(1:nasusc2)
  if (itc == 1) write(*,'("************************************************************",/)')
! All done!
  end subroutine check_atomic_options
