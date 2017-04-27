  subroutine set_spinrot_options(itc,my_rank)
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
  integer(kind=i4b) :: ig, istart, iend
  real(kind=r8b)    :: ulen


! No spin rotations will be done
  if (.not.lrot) then
    if (my_rank == 0) then 
      if (itc == 1) then
        write(*,'("************************************************************")')
        write(*,'(" No spin rotations")')
        write(*,'("************************************************************"/)')
      end if
    end if ! my_rank 
    return
  end if

! Spin rotations options

! rigid rotation or atom dependent
  call find_keyinfile('ispinrot',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) ispinrot
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_spinrot_options: WARNING - key ispinrot not found!")')
    end if ! my_rank
    ispinrot = 0
  end if
  if (ispinrot /= 0 .and. ispinrot /= 1)  stop 'set_spinrot_options: unknown spin rotation option'

! new global spin direction
  call find_keyinfile('urot',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) urot(1:3)
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_spinrot_options: WARNING - key urot not found!")')
    end if ! my_rank
    urot = (/ 0.d0, 0.d0, 1.d0 /)
  end if
  ulen = sqrt(dot_product(urot,urot))
  if (ulen < 1.d-6) stop 'set_spinrot_options: urot has length zero!'
  urot = urot/ulen

! mix factor for new scf orientations from torques
  call find_keyinfile('dirmix2',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) dirmix2
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_spinrot_options: WARNING - key dirmix2 not found!")')
    end if ! my_rank 
    dirmix2 = 0.d0
  end if

! mix factor for new scf orientations from output orientations
  call find_keyinfile('dirmix3',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) dirmix3
  else
    if (my_rank == 0) then 
      if (lwarn .and. itc == 1) write(*,'("set_spinrot_options: WARNING - key dirmix3 not found!")')
    end if ! my_rank
    dirmix3 = 0.5d0
  end if

! input GF is in collinear ASA form
  call find_keyinfile('lgrefsph',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lgrefsph
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_spinrot_options: WARNING - key lgrefsph not found!")')
    end if ! my_rank
    lgrefsph = .false.
  end if

! ----------------------------------------------------------------------
! Print options summary
  if (my_rank == 0) then 
    if (itc == 1) then
      write(*,'("************************************************************")')
      if (ispinrot == 0) write(*,'(" Global spin rotation: urot=",3f8.4)') urot
      if (ispinrot == 1) write(*,'(" Local spin rotations: check magdir.dat!")')
      if (ispinrot == 1) write(*,'(" Mix from torques with dirmix2=",f8.4)') dirmix2
      if (ispinrot == 1) write(*,'(" Mix from outputs with dirmix3=",f8.4)') dirmix3
      if (lgrefsph) then
        write(*,'(" Input GF assumed to be collinear ASA  --  fast spin rotations")')
      else
        write(*,'(" Input GF assumed to be non-collinear FP  --  slow spin rotations")')
      end if
      write(*,'("************************************************************"/)')
    end if
  end if ! my_rank 
! All done!
  end subroutine set_spinrot_options
