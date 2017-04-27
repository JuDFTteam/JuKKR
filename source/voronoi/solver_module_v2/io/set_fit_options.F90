  subroutine set_fit_options(itc,my_rank)
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
  real(kind=r8b)    :: re, im


! No rational fit of GF
  if (.not.lfit) then
    if (my_rank == 0) then
      if (itc == 1) then
        write(*,'("************************************************************")')
        write(*,'(" No fitting of GF")')
        write(*,'("************************************************************"/)')
      end if
    end if ! my_rank 
    return
  end if

! Rational fit options

! type of fit
  call find_keyinfile('ifit',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) ifit
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_fit_options: WARNING - key ifit not found!")')
    end if ! my_rank
    ifit = 2
  end if
  if (ifit /= 1 .and. ifit /= 2) stop 'set_fit_options: unknown fit/interpolation option'

! degree of numerator of rational fit function
  call find_keyinfile('numd',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) numd
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_fit_options: WARNING - key numd not found!")')
    end if ! my_rank
    numd = 19
  end if

! degree of denominator of rational fit function
  call find_keyinfile('dend',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) dend
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_fit_options: WARNING - key dend not found!")')
    end if ! my_rank
    dend = 20
  end if

! energy zero for fit polynomials
  call find_keyinfile('eshift',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) re, im
    eshift = cmplx(re,im)
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_fit_options: WARNING - key eshift not found!")')
    end if ! my_rank
    eshift = 0.d0
  end if

! whether to shift real part of GF
  call find_keyinfile('lregf',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lregf
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_fit_options: WARNING - key lregf not found!")')
    end if ! my_rank
    lregf = .false.
  end if

! by how much to shift real part of GF
  call find_keyinfile('fudge',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) fudge
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_fit_options: WARNING - key fudge not found!")')
    end if ! my_rank
    fudge = 0.d0
  end if

! Print options summary
  if (my_rank == 0) then
    if (itc == 1) then
      write(*,'("************************************************************")')
      write(*,'(" Fit to energy dependence of GF")')
      if (ifit == 1) write(*,'(" + GF interpolated with ndeg=",i4)') numd
      if (ifit == 2) write(*,'(" + GF fit to rational function with numd, dend=",2i4)') numd, dend
      if (abs(eshift) > 1.d-8) write(*,'(" + Energy zero in polynomials shifted to",2f12.6)') eshift
!     shift of Re GF
      if (lregf) write(*,'(" + Shift of Re GF, fudge=",f10.6)') fudge
      write(*,'("************************************************************"/)')
    end if
  end if ! my_rank
! All done!
  end subroutine set_fit_options
