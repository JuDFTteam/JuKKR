  subroutine set_susc_options(itc,my_rank)
! global options read in from input file
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
  integer(kind=i4b) :: ig, istart, iend, i, j
  real(kind=r8b)    :: ere, eim


! Not a susceptibility calculation
  if (.not.lsusc) then
!   Set overrides
    lkha = .false.; lkxc = .false.
    if (my_rank == 0) then
      if (itc == 1) write(*,'("************************************************************")')
      if (itc == 1) write(*,'(" No susceptibility calculation")')
      if (itc == 1) write(*,'("************************************************************"/)')
    end if ! my_rank
    return
  end if

! Susceptibility options

! static or dynamic susceptibility
  call find_keyinfile('ldynsusc',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) ldynsusc
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_susc_options: WARNING - key ldynsusc not found!")')
      ldynsusc = .false.
    end if ! my_rank
  end if

! how many frequencies
  call find_keyinfile('nomega',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) nomega
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_susc_options: WARNING - key nomega not found!")')
    end if ! my_rank
    nomega = 201
  end if

! first frequency
  call find_keyinfile('omegamin',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) ere, eim
    omegamin = cmplx(ere,eim)
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_susc_options: WARNING - key omegamin not found!")')
    end if ! my_rank 
    omegamin = (-1.d-3,0.d0)
  end if

! last frequency
  call find_keyinfile('omegamax',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) ere, eim
    omegamax = cmplx(ere,eim)
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_susc_options: WARNING - key omegamax not found!")')
    end if ! my_rank
    omegamax = (+1.d-3,0.d0)
  end if

! frequency increment
  call find_keyinfile('domega',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) domega
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_susc_options: WARNING - key domega not found!")')
    end if ! my_rank
    domega = 0.d0
  end if

! cutoff for spherical harmonic expansion
  call find_keyinfile('nlmax0',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) nlmax0
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_susc_options: WARNING - key nlmax0 not found!")')
    end if ! my_rank
    nlmax0 = 2*nlmax
  end if
  lmmax0 = (nlmax0 + 1)**2               ! lm components for susceptibility

! include Hartree kernel
  call find_keyinfile('lkha',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lkha
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_susc_options: WARNING - key lkha not found!")')
    end if ! my_rank 
    lkha = .false.
  end if

! include xc kernel
  call find_keyinfile('lkxc',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lkxc
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_susc_options: WARNING - key lkxc not found!")')
    end if ! my_rank
    lkxc = .true.
  end if

! output representation for susceptibility
  call find_keyinfile('lcartesian',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lcartesian
  else 
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_susc_options: WARNING - key lcartesian not found!")')
    end if ! my_rank
    lcartesian = .true.
  end if

! include analytic integral
  call find_keyinfile('lanalytic',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lanalytic
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_susc_options: WARNING - key lanalytic not found!")')
    end if ! my_rank 
    lanalytic = .true.
  end if

! include nonanalytic integral
  call find_keyinfile('lnonanalytic',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lnonanalytic
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_susc_options: WARNING - key lnonanalytic not found!")')
    end if ! my_rank
    lnonanalytic = .true.
  end if

! solve Dyson equation for susceptibility
  call find_keyinfile('lenhanced',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lenhanced
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_susc_options: WARNING - key lenhanced not found!")')
    end if ! my_rank
    lenhanced = .true.
  end if

! iterations to try to get Goldstone mode
  call find_keyinfile('itermax',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) itermax
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_susc_options: WARNING - key itermax not found!")')
    end if ! my_rank
    itermax = 0
  end if

! mix factor for Goldstone mode
  call find_keyinfile('lambdamix',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lambdamix
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_susc_options: WARNING - key lambdamix not found!")')
    end if ! my_rank 
    lambdamix = 0.d0
  end if

! use spin sum rule
  call find_keyinfile('lsumrule',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lsumrule
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_susc_options: WARNING - key lsumrule not found!")')
    end if ! my_rank
    lsumrule = .true.
  end if

! use noncollinear spin sum rule
  call find_keyinfile('lnewsumrule',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lnewsumrule
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_susc_options: WARNING - key lnewsumrule not found!")')
    end if ! my_rank
    lnewsumrule = .true.
  end if

! use onsite gf in susc
  call find_keyinfile('lsusconsite',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lsusconsite
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_susc_options: WARNING - key lsusconsite not found!")')
    end if ! my_rank
    lsusconsite = .true.
  end if

! use structural gf in susc
  call find_keyinfile('lsuscstruct',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lsuscstruct
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_susc_options: WARNING - key lsuscstruct not found!")')
    end if ! my_rank
    lsuscstruct = .true.
  end if

! spin-current spin correlation function (added by Sascha)
  call find_keyinfile('lcurrcorr',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lcurrcorr
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_susc_options: WARNING - key lcurrcorr not found!")')
    end if ! my_rank
    lcurrcorr = .false.
  end if

! Interpolation of spin-current spin correlation function (added by Sascha)
  call find_keyinfile('lcurrcorrint',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lcurrcorrint
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_susc_options: WARNING - key lcurrcorrint not found!")')
    end if ! my_rank
    lcurrcorrint = .false.
  end if

! Divergence of spin-current function (added by Sascha)
  call find_keyinfile('lcurrcorrdiv',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lcurrcorrdiv
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_susc_options: WARNING - key lcurrcorrdiv not found!")')
    end if ! my_rank
    lcurrcorrdiv = .false.
  end if

! Scalar relativistic correction to the continuity equation
  call find_keyinfile('lscalarcorr',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lscalarcorr
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_susc_options: WARNING - key lscalarcorr not found!")')
    end if ! my_rank
    lscalarcorr = .false.
  end if

! Defaults
  isusc(1:nasusc) = 0; ikha(1:nasusc) = 0; ikxc(1:nasusc) = 0
  nasusc2 = 0; iasusc2 = 0

! Search for atom info lines
! first pass: which atoms included in susceptibility calculation
  do iline=1,nlines
    call find_keyinline('ig',nchars,nlines,inputfile,iline,ipos,found)
    if (found) then
      read(inputfile(iline)(ipos:nchars),*) ig
      if (ig > ngroup) cycle
!     which atoms are in group ig
      istart = sum(igroup(1:ig-1)) + 1
      iend   = sum(igroup(1:ig))
!     susceptibility for this atom group
      call find_keyinline('isusc',nchars,nlines,inputfile,iline,ipos,found)
      if (found) then
        read(inputfile(iline)(ipos:nchars),*) i
        if (i == 0) cycle
        if (i /= 1 .and. i /= 2 .and. i /= 3) stop 'set_susc_options: unknown isusc option'
        isusc(istart:iend) = i
!       update list of susc atoms
        do j=sum(igroup(1:ig-1))+1,sum(igroup(1:ig))
          iasusc2(j) = j
        end do
        nasusc2 = nasusc2 + igroup(ig)
      end if
    end if
  end do
! second pass: atomic susceptibility options
  do iline=1,nlines
    call find_keyinline('ig',nchars,nlines,inputfile,iline,ipos,found)
    if (found) then
      read(inputfile(iline)(ipos:nchars),*) ig
      if (ig > ngroup) cycle
!     which atoms are in group ig
      istart = sum(igroup(1:ig-1)) + 1
      iend   = sum(igroup(1:ig))
!     Hartree kernel
      if (lkha .and. isusc(istart) /= 0) then
        call find_keyinline('ikha',nchars,nlines,inputfile,iline,ipos,found)
        if (found) then
          read(inputfile(iline)(ipos:nchars),*) i
          if (i /= 0 .and. i /= 1) stop 'set_susc_options: unknown ikha option'
          ikha(istart:iend) = i
        end if
      end if
!     xc kernel
      if (lkxc .and. isusc(istart) /= 0) then
        call find_keyinline('ikxc',nchars,nlines,inputfile,iline,ipos,found)
        if (found) then
          read(inputfile(iline)(ipos:nchars),*) i
          if (i /= 0 .and. i /= 1 .and. i /= 2 .and. i /= 3) stop 'set_susc_options: unknown ikxc option'
          ikxc(istart:iend) = i
        end if
      end if
    end if
  end do

! Print options summary
  if (my_rank == 0) then
    if (itc == 1) then
      write(*,'("************************************************************")')
      write(*,'(" Susceptibility    ==>  check atomic options")')
!     angular momentum cutoff
      write(*,'(" + lmax for susceptibility is",i4)') nlmax0
!     static or dynamic, which integrals
      if (ldynsusc) then
        write(*,'(" + Dynamic calculation:")')
        write(*,'(" + nomega, omegamin, omegamax, domega=",i6,5f10.6)') nomega, omegamin, omegamax, domega
        write(*,'(" + Analytic integral=",l1,", non-analytic integral=",l1)') lanalytic, lnonanalytic
      else
        write(*,'(" + Static calculation")')
      end if
!     Kohn-Sham only or enhanced
      if (lenhanced) then
        write(*,'(" + Enhanced susceptibility")')
!       spin sum rule
        if (lsumrule) then
          if (lnewsumrule) then
            write(*,'(" + Spin sum rule for noncollinear case")')
          else
            write(*,'(" + Spin sum rule for collinear case")')
          end if
        else
          write(*,'(" + Spin sum rule not used ==> DANGER !!!")')
        end if
      else
        write(*,'(" + Kohn-Sham susceptibility only")')
      end if
!     output format
      if (lcartesian) then
        write(*,'(" + Susceptibility in cartesian components")')
      else
        write(*,'(" + Susceptibility in spin components")')
      end if
!     shake your Hartree
      if (lkha) then
        write(*,'(" Hartree kernel    ==>  check atomic options")')
        if (.not.lpositions) stop 'positions must be provided in separate file'
      else
        write(*,'(" No Hartree kernel computed")')
      end if
!     to xc or not to xc, that is the question
      if (lkxc) then
        write(*,'(" xc kernel         ==>  check atomic options")')
      else
        write(*,'(" No xc kernel computed")')
      end if
!     Spin-currents
      if (lcurrcorr) then
        write(*,'(" + Spin-current spin-density correlation function")')
        if(lcurrcorrint) then
          write(*,'("   Interpolation on grid of size",i3)') n_int
        end if
        if(lcurrcorrdiv) then
          write(*,'("   Divergence calculated (be carful because of numerical problems)")')
        end if
      end if
!     Info on atoms selected for susc calculation
      write(*,'(/" Number of atoms for susceptibility calculation:",i4)') nasusc2
      write(*,'(" iasusc2=",1000i4)') iasusc2(1:nasusc2)
      write(*,'(/" Susceptibility configuration for each atom group:")')
      do ig=1,ngroup
        istart = sum(igroup(1:ig-1)) + 1
        iend   = sum(igroup(1:ig))
        if (isusc(istart) /= 0) then
          write(*,'(" ig, istart, iend=",3i6,"  iasusc=",1000i4)') ig, istart, iend, iasusc(istart:iend)
!         Which part of the susceptibility to compute
          if (isusc(istart) == 1) write(*,'(" + Transverse susceptibility")')
          if (isusc(istart) == 2) write(*,'(" + Longitudinal susceptibility")')
          if (isusc(istart) == 3) write(*,'(" + Full susceptibility")')
!         Hartree kernel option
          if (lkha) then
            if (ikha(istart) == 0)  write(*,'(" + No Hartree kernel")')
            if (ikha(istart) == 1)  write(*,'(" + Hartree kernel")')
          end if
!         xc kernel options
          if (lkxc) then
            if (ikxc(istart) == 0)  write(*,'(" + No xc kernel")')
            if (ikxc(istart) == 1)  write(*,'(" + Transverse xc kernel")')
            if (ikxc(istart) == 2)  write(*,'(" + Longitudinal xc kernel")')
            if (ikxc(istart) == 3)  write(*,'(" + Full xc kernel")')
          end if
        end if
      end do
      write(*,'("************************************************************"/)')
    end if
  end if ! my_rank 
! All done!
  end subroutine set_susc_options
