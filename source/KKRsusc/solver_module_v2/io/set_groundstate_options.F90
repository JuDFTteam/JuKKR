  subroutine set_groundstate_options(itc,my_rank)
! groundstate options read in from input file
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
  integer(kind=i4b) :: ig, istart, iend, ntot, i


! Groundstate options

! diagonalize atomic density matrix
  call find_keyinfile('lrhodiag',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lrhodiag
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_groundstate_options: WARNING - key lrhodiag not found!")')
    end if ! my_rank
    lrhodiag = .false.
  end if

! exclude onsite GF from GS properties
  call find_keyinfile('lgsonsite',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lgsonsite
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_groundstate_options: WARNING - key lgsonsite not found!")')
    end if ! my_rank
    lgsonsite = .true.
  end if

! exclude structural GF from GS properties
  call find_keyinfile('lgsstruct',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lgsstruct
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_groundstate_options: WARNING - key lgsstruct not found!")')
    end if ! my_rank
    lgsstruct = .true.
  end if

! update charge density
  call find_keyinfile('lnewdensity',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lnewdensity
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_groundstate_options: WARNING - key lnewdensity not found!")')
    end if ! my_rank
    lnewdensity = .true.
  end if

! read the number of energy points for SCF integration from the file
  call find_keyinfile('lscfmesh',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lscfmesh
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_groundstate_options: WARNING - key lscfmesh not found!")')  
    end if ! my_rank
!   mesh must be provided separately when using fitted GF
    if (lfit) then
      lscfmesh = .true.
    else
      lscfmesh = .false.
    end if
  end if

! read number of energy points from file?
  if (lscfmesh) then
    open(file='emesh.scf',unit=iofile,status='old')
    read(iofile,*) nescf
    close(iofile)
  else
    nescf = nesusc
  end if

! replace structural GF with free space GF
  call find_keyinfile('lfreegf',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lfreegf
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_groundstate_options: WARNING - key lfreegf not found!")')
    end if ! my_rank
    lfreegf = .false.
  end if

! whether to include the Fermi energy in the definition of the band energy
  call find_keyinfile('lebandnoef',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lebandnoef
  else
    if (my_rank == 0) then 
      if (lwarn .and. itc == 1) write(*,'("set_groundstate_options: WARNING - key lebandnoef not found!")')
    end if ! my_rank
    lebandnoef = .false.
  end if

! read atomic positions from file
  call find_keyinfile('lpositions',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lpositions
  else
    if (my_rank == 0) then 
      if (lwarn .and. itc == 1) write(*,'("set_groundstate_options: WARNING - key lpositions not found!")')
    end if ! my_rank
    lpositions = .false.
  end if

! kill spin density in selected atoms
  call find_keyinfile('lnobxc',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lnobxc
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_groundstate_options: WARNING - key lnobxc not found!")')
    end if ! my_rank 
    lnobxc = .false.
    inobxc(1:nasusc) = 0
  end if

! Search for atom info lines
  if (lnobxc) then
    ntot = 0
    do iline=1,nlines
      call find_keyinline('ig',nchars,nlines,inputfile,iline,ipos,found)
      if (found) then
        read(inputfile(iline)(ipos:nchars),*) ig
        if (ig > ngroup) cycle
        istart = sum(igroup(1:ig-1)) + 1
        iend   = sum(igroup(1:ig))
!       which atoms are in group ig
        call find_keyinline('inobxc',nchars,nlines,inputfile,iline,ipos,found)
        if (found) then
          read(inputfile(iline)(ipos:nchars),*) i
          if (i /= 0) then
            if (i /= 1) stop 'set_groundstate_options: unknown inobxc option'
            inobxc(istart:iend) = i
            ntot = ntot + igroup(ig)
          else
            inobxc(istart:iend) = 0
          end if
        end if
      end if
    end do
  end if

! calculate groundstate current
  call find_keyinfile('lcurrent',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lcurrent
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_groundstate_options: WARNING - key lcurrent not found!")')
    end if ! my_rank 
    lcurrent = .false.
  end if

! add soc contribution to the current
  call find_keyinfile('lcurrentsoc',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lcurrentsoc
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_groundstate_options: WARNING - key lcurrentsoc not found!")')
    end if ! my_rank 
    lcurrentsoc = .false.
  end if

! add zeeman contribution to the current
  call find_keyinfile('lcurrentzeeman',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lcurrentzeeman
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_groundstate_options: WARNING - key lcurrentzeeman not found!")')
    end if ! my_rank 
    lcurrentzeeman = .false.
  end if

! output currents to files
  call find_keyinfile('lcurrentoutput',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lcurrentzeeman
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_groundstate_options: WARNING - key lcurrentoutput not found!")')
    end if ! my_rank 
    lcurrentoutput = .false.
  end if

! Interpolation of groundstate current
  call find_keyinfile('lcurrentint',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lcurrentint
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_groundstate_options: WARNING - key lcurrentint not found!")')
    end if ! my_rank 
    lcurrentint = .false.
  end if

! Current induced magnetic fields
  call find_keyinfile('lcurrentbfield',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lcurrentbfield
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_groundstate_options: WARNING - key lcurrentbfield not found!")')
    end if ! my_rank 
    lcurrentbfield = .false.
  end if

! Interpolation grid size
! used for plotting of the interpolated currents
  call find_keyinfile('nint',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) n_int
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_groundstate_options: WARNING - key n_int not found!")')
    end if ! my_rank 
    n_int = 25
  end if

! Print options summary
  if (my_rank == 0) then
    if (itc == 1) then
      write(*,'("************************************************************")')
      write(*,'(" Groundstate options")')
      if (lrhodiag)  write(*,'(" + Atomic density matrix diagonalized")')
      if (.not.lgsonsite) write(*,'(" + Onsite GF not included in groundstate properties")')
      if (.not.lgsstruct) write(*,'(" + Structural GF not included in groundstate properties")')
      if (lnewdensity)    write(*,'(" + Charge density updated")')
      if (lscfmesh)       write(*,'(" + SCF mesh read from file")')
      if (lfreegf)        write(*,'(" + Structural GF replaced by free space GF")')
      if (lebandnoef)     write(*,'(" + Band energy does not include Fermi energy")')
      if (lpositions)     write(*,'(" + Atomic positions read from file")')
      if (lcurrent)       write(*,'(" + Groundstate currents (paramagnetic)")')
      if (lcurrentint)    write(*,'(" + Interpolation of groundstate currents with n_int= ",i3)') n_int
      if (lcurrentsoc)    write(*,'(" + SOC contribution to groundstate currents ")')
      if (lcurrentzeeman) write(*,'(" + Zeeman contribution to groundstate charge currents ")')
!     atom info
      if (lnobxc) then
        write(*,'(" + Kill spin density in selected atoms")')
        write(*,'(/" Number of atoms to kill:",i4)') ntot
        write(*,'(/" Kill configuration for each atom group:")')
        do ig=1,ngroup
          istart = sum(igroup(1:ig-1)) + 1
          iend   = sum(igroup(1:ig))
          if (inobxc(istart) /= 0) then
            write(*,'(" ig, istart, iend=",3i6,"  iasusc=",1000i4)') ig, istart, iend, iasusc(istart:iend)
          end if
        end do
      end if
      write(*,'("************************************************************",/)')
    end if
  end if ! my_rank 
! All done!
  end subroutine set_groundstate_options
