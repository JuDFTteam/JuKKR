  subroutine set_basis_options(itc,my_rank)
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
  integer(kind=i4b) :: ig, istart, iend, is, il, ib, iw(0:nlmax)
  real(kind=r8b)    :: ere1, eim1, ere2, eim2
  logical           :: lbasisdefault


! Basis construction options

! basis type
  call find_keyinfile('ibasis',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) ibasis
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_basis_options: WARNING - key ibasis not found!")')
    end if ! my_rank
    ibasis = 2
  end if
  if (ibasis /=1 .and. ibasis /= 2 .and. ibasis /= 3) stop 'set_basis_options: unknown basis option'

! construction method
  call find_keyinfile('ibasismethod',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) ibasismethod
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_basis_options: WARNING - key ibasismethod not found!")')
    end if ! my_rank
    ibasismethod = 2
  end if
  if (ibasismethod /= 0 .and. ibasismethod /= 1 .and. ibasismethod /= 2) stop 'set_basis_options: unknown basis orthogonalization option'

! significance threshold
  call find_keyinfile('basistol',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) basistol
  else
    if (my_rank == 0) then 
      if (lwarn .and. itc == 1) write(*,'("set_basis_options: WARNING - key basistol not found!")')
    end if ! my_rank
    basistol = 1.d-3
  end if

! keep basis fixed or not
  call find_keyinfile('lbasisfreeze',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lbasisfreeze
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_basis_options: WARNING - key lbasisfreeze not found!")')
    end if ! my_rank
    lbasisfreeze = .false.
  end if

! project t-matrix in basis
  call find_keyinfile('ltmatproj',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) ltmatproj
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_basis_options: WARNING - key ltmatproj not found!")')
    end if ! my_rank
    ltmatproj = .false.
  end if

!
! use default setup
  call find_keyinfile('lbasisdefault',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lbasisdefault
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_basis_options: WARNING - key lbasisdefault not found!")')
    end if ! my_rank
    lbasisdefault = .false.
  end if

! Use defaults or read everything from input file

! ***********************
  if (lbasisdefault) then
! ***********************

!   default construction method
    ibasis = 2; ibasismethod = 2; basistol = 1.d-3
!   default projection energies
    issusc(1:nasusc) = nsmax
    iwsusc(0:nlmax,1:nsmax,1:nasusc) = 4
    ewsusc(1,0:nlmax,1:nsmax,1:nasusc) = (0.3d0,0.d0)
    ewsusc(2,0:nlmax,1:nsmax,1:nasusc) = (0.5d0,0.d0)
    ewsusc(3,0:nlmax,1:nsmax,1:nasusc) = (0.7d0,0.d0)
    ewsusc(4,0:nlmax,1:nsmax,1:nasusc) = (0.9d0,0.d0)

! ****
  else
! ****

!   search for atom info lines
!   first pass: find out basis structure for each atom group
!   +++++++++++++++++
    do iline=1,nlines
!   +++++++++++++++++
!     search for atom info lines
      call find_keyinline('ig',nchars,nlines,inputfile,iline,ipos,found)
      if (found) then
        read(inputfile(iline)(ipos:nchars),*) ig
        istart = sum(igroup(1:ig-1)) + 1
        iend   = sum(igroup(1:ig))
!       which atoms are magnetic (is == 2) or not (is == 1)
        call find_keyinline('nspin',nchars,nlines,inputfile,iline,ipos,found)
        if (found) then
          read(inputfile(iline)(ipos:nchars),*) is
          if (is > nsmax) stop 'set_basis_options: is > nsmax'
          issusc(istart:iend) = is
        end if
!       info on how many basis functions for each l channel
        call find_keyinline('iwsusc',nchars,nlines,inputfile,iline,ipos,found)
        if (found) then
          read(inputfile(iline)(ipos:nchars),*) iw(0:nlmax)
          do il=0,nlmax
            iwsusc(il,1:nsmax,istart:iend) = iw(il)
          end do
        end if
      end if
!   ++++++
    end do
!   ++++++
!   second pass: read projection energies
!   +++++++++++++++++
    do iline=1,nlines
!   +++++++++++++++++
!     search for atom info lines
      call find_keyinline('ig',nchars,nlines,inputfile,iline,ipos,found)
      if (found) then
        read(inputfile(iline)(ipos:nchars),*) ig
        if (ig > ngroup) cycle
        istart = sum(igroup(1:ig-1)) + 1
        iend   = sum(igroup(1:ig))
!       find l-channel
        call find_keyinline('il',nchars,nlines,inputfile,iline,ipos,found)
        if (found) then
          read(inputfile(iline)(ipos:nchars),*) il
          if (il > nlmax) stop 'set_basis_options: il > nlmax'
!         find projection energy
          call find_keyinline('ib',nchars,nlines,inputfile,iline,ipos,found)
          if (found) then
            read(inputfile(iline)(ipos:nchars),*) ib
            if (ib > iwsusc(il,1,istart)) stop 'set_basis_options: ib > nbmax'
            call find_keyinline('ewsusc',nchars,nlines,inputfile,iline,ipos,found)
            if (.not.found) stop 'set_basis_options: projection energy not found'
!           one or two projection energies
            if (issusc(istart) == 1) then
              read(inputfile(iline)(ipos:nchars),*) ere1, eim1
              ewsusc(ib,il,1,istart:iend) = cmplx(ere1,eim1)
            else
              read(inputfile(iline)(ipos:nchars),*) ere1, eim1, ere2, eim2
              ewsusc(ib,il,1,istart:iend) = cmplx(ere1,eim1)
              ewsusc(ib,il,2,istart:iend) = cmplx(ere2,eim2)
            end if
          end if
        end if
      end if
!   ++++++
    end do
!   ++++++

! ******
  end if
! ******

! Check basis options
  if (my_rank == 0) then
    if (itc == 1) then
      write(*,'("************************************************************")')
      write(*,'(" Basis construction")')
      if (lbasisdefault) then
        write(*,'(" + Using default options")')
      else
        write(*,'(" + Using input options  =>  check atomic options")')
       end if
      if (ltmatproj) then
        write(*,'(" + Starting t-matrix computed in projection basis")')
      else
        write(*,'(" + Starting t-matrix read in from host KKR program")')
      end if

      if (lbasisfreeze) then
        write(*,'(" + Basis functions read in from file")')
      else
        if (ibasis == 1) write(*,'(" + Basis functions for each l-channel and spin")')
        if (ibasis == 2) write(*,'(" + Basis functions for each l-channel, common to both spins")')
        if (ibasis == 3) write(*,'(" + Basis functions common for all l-channels and spins")')
      end if
      if (ibasismethod == 0) write(*,'(" + No basis reduction performed, for testing only!")')
      if (ibasismethod == 1) write(*,'(" + Gram-Schmidt orthogonalization, tol=",f10.6)') basistol
      if (ibasismethod == 2) write(*,'(" + Diagonalization of overlap matrix, tol=",f10.6)') basistol
      if (.not.lbasisdefault) then
        write(*,'(/"  Basis configuration for each atom group:")')
        do ig=1,ngroup
          istart = sum(igroup(1:ig-1)) + 1
          iend   = sum(igroup(1:ig))
          write(*,'(" ig, istart, iend=",3i6,"  iasusc=",1000i4)') ig, istart, iend, iasusc(istart:iend)
          write(*,'(" + issusc=",i2)') issusc(istart)
          do il=0,nlmax
            write(*,'(" + il=",i2,"  iwsusc=",2i4)') il, iwsusc(il,1:nsmax,istart)
            do ib=1,maxval(iwsusc(il,1:nsmax,istart))
              write(*,'("     ib=",i4,"  ewsusc=",4f8.4)') ib, ewsusc(ib,il,1:nsmax,istart)
            end do
          end do
        end do
      end if
      write(*,'("************************************************************"/)')
    end if
  end if ! my_rank
! All done!
  end subroutine set_basis_options
