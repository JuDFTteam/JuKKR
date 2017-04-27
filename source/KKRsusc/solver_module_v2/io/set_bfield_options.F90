  subroutine set_bfield_options(itc,my_rank)
! bfield options read in from input file
  use global

  implicit none


! Current iteration
  integer(kind=i4b), intent(in) :: itc, my_rank
! ----------------------------------------------------------------------
! line number, column number
  integer(kind=i4b) :: iline, ipos
! was key found?
  logical           :: found
! ----------------------------------------------------------------------
  integer(kind=i4b) :: ig, istart, iend, i, ntot, ibfield0
  logical           :: lbfield0
  real(kind=r8b)    :: blen0, bdir0(1:3)


! No magnetic field applied
  if (.not.lbfield) then 
    if (my_rank == 0) then 
      if (itc == 1) then
        write(*,'("************************************************************")')
        write(*,'(" No external B field")')
        write(*,'("************************************************************"/)')
      end if
    end if ! my_rank
    return
  end if

! Magnetic field options

! use constraining fields
  call find_keyinfile('lbconst',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lbconst
    if (lbconst .and. .not. ljij) stop 'set_bfield_options: constraining fields need Jij!'
    if (lbconst .and. .not. lrot) stop 'set_bfield_options: constraining fields need spin rotations!'
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_bfield_options: WARNING - key lbconst not found!")')
    end if ! my_rank
    lbconst = .false.
  end if

! set defaults
  call find_keyinfile('lbfield0',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lbfield0
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_bfield_options: WARNING - key lbfield0 not found!")')
    end if ! my_rank
    lbfield0 = .false.
  end if

! ********************
  if (lbfield0) then
! ********************

!   Same magnetic field for all atoms

!   form of Zeeman coupling
    call find_keyinfile('ibfield0',nchars,nlines,inputfile,iline,ipos,found)
    if (found) then
      read(inputfile(iline)(ipos:nchars),*) ibfield0
    else
      if (my_rank == 0) then
        if (lwarn .and. itc == 1) write(*,'("set_bfield_options: WARNING - key ibfield0 not found!")')
      end if ! my_rank 
      ibfield0 = 1
    end if

!   magnetic field intensity
    call find_keyinfile('blen0',nchars,nlines,inputfile,iline,ipos,found)
    if (.not.found) stop 'set_bfield_options: key blen0 not found!'
    read(inputfile(iline)(ipos:nchars),*) blen0

!   magnetic field direction
    call find_keyinfile('bdir0',nchars,nlines,inputfile,iline,ipos,found)
    if (found) then
      read(inputfile(iline)(ipos:nchars),*) bdir0(1:3)
      blen0 = sqrt(dot_product(bdir0,bdir0))
      if (blen0 < 1.d-6) stop 'set_bfield_options: bdir0 has length zero'
      bdir0 = bdir0/blen0
    else
      if (my_rank == 0) then
        if (lwarn .and. itc == 1) write(*,'("set_bfield_options: WARNING - key bdir0 not found!")')
      end if ! my_rank
      bdir0 = (/ 0.d0, 0.d0, 1.d0 /)
    end if

!   all atoms set to the same thing
    ntot = nasusc
    ibfield(1:nasusc) = ibfield0; blen(1:nasusc) = blen0
    bdir(:,1:nasusc)  = spread(bdir0,2,nasusc)

! ****
  else
! ****

!   Read atom dependent magnetic fields

!   Defaults
    ntot = 0
    ibfield(1:nasusc) = 0; blen(1:nasusc) = 0.d0
    bdir(:,1:nasusc) = spread((/ 0.d0, 0.d0, 1.d0 /),2,nasusc)

!   Search for atom info lines
!   first pass: which atoms will have applied B field
    do iline=1,nlines
      call find_keyinline('ig',nchars,nlines,inputfile,iline,ipos,found)
      if (found) then
        read(inputfile(iline)(ipos:nchars),*) ig
        if (ig > ngroup) cycle
        istart = sum(igroup(1:ig-1)) + 1
        iend   = sum(igroup(1:ig))
!       local form of Zeeman coupling
        call find_keyinline('ibfield',nchars,nlines,inputfile,iline,ipos,found)
        if (found) then
          read(inputfile(iline)(ipos:nchars),*) i
          if (i /= 0) then
            if (i /= 1 .and. i /= 2 .and. i /= 3) stop 'set_bfield_options: unknown ibfield option'
            ibfield(istart:iend) = i
            ntot = ntot + igroup(ig)
          end if
        end if
      end if
    end do
!   second pass: atomic B field options
    do iline=1,nlines
      call find_keyinline('ig',nchars,nlines,inputfile,iline,ipos,found)
      if (found) then
        read(inputfile(iline)(ipos:nchars),*) ig
        if (ig > ngroup) cycle
        istart = sum(igroup(1:ig-1)) + 1
        iend   = sum(igroup(1:ig))
!       magnetic field intensity
        call find_keyinline('blen',nchars,nlines,inputfile,iline,ipos,found)
        if (found .and. ibfield(istart) /= 0) then
          read(inputfile(iline)(ipos:nchars),*) blen0
          blen(istart:iend) = blen0
        end if
!       magnetic field direction
        call find_keyinline('bdir',nchars,nlines,inputfile,iline,ipos,found)
        if (found .and. ibfield(istart) /= 0) then
          read(inputfile(iline)(ipos:nchars),*) bdir0(1:3)
          blen0 = sqrt(dot_product(bdir0,bdir0))
          if (blen0 < 1.d-6) stop 'set_bfield_options: bdir has length zero'
          bdir0 = bdir0/blen0
          bdir(:,istart:iend) = spread(bdir0,2,igroup(ig))
        end if
      end if
    end do

! ******
  end if
! ******

! Print options summary
  if (my_rank == 0) then 
    if (itc == 1) then
      write(*,'("************************************************************")')
!     same for all atoms
      if (lbfield0) then
        write(*,'(" Uniform external B field")')
        if (ibfield0 == 1) write(*,'(" B field in B.S form, B=",3es12.4)') blen(1)*bdir(:,1)
        if (ibfield0 == 2) write(*,'(" B field in B.L form, B=",3es12.4)') blen(1)*bdir(:,1)
        if (ibfield0 == 3) write(*,'(" B field in B.(L+S) form, B=",3es12.4)') blen(1)*bdir(:,1)
!       if (ibfield0 == 4) write(*,'(" Constraining field, magdir=",4es12.4)') blen(1), bdir(:,1)
        return
      end if
!     atom-specific magnetic fields
      write(*,'(" Atom-dependent external B field  ==>  check atomic options")')
!     Atom info
      write(*,'(/" Number of atoms with B field:",i4)') ntot
      write(*,'(/" B field configuration for each atom group:")')
      do ig=1,ngroup
        istart = sum(igroup(1:ig-1)) + 1
        iend   = sum(igroup(1:ig))
        if (ibfield(istart) /= 0) then
          write(*,'(" ig, istart, iend=",3i6,"  iasusc=",1000i4)') ig, istart, iend, iasusc(istart:iend)
          if (ibfield(istart) == 1) write(*,'(" + B field in B.S form, B=",3es12.4)') blen(istart)*bdir(:,istart)
          if (ibfield(istart) == 2) write(*,'(" + B field in B.L form, B=",3es12.4)') blen(istart)*bdir(:,istart)
          if (ibfield(istart) == 3) write(*,'(" + B field in B.(L+S) form, B=",3es12.4)') blen(istart)*bdir(:,istart)
!         if (ibfield(istart) == 4) write(*,'(" + Constraining field, magdir=",4es12.4)') blen(istart), bdir(:,istart)
        end if
      end do
      if (lbconst) write(*,'(/," Constraining B fields  ==>  check Jij options")')
      write(*,'("************************************************************"/)')
    end if
  end if ! my_rank
! All done!
  end subroutine set_bfield_options
