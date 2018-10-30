  subroutine load_input(filename,nchars,nlines,my_rank)
! reads input file to memory
  use global, only: i4b, iomain, iofile, inputfile

  implicit none


! name of input file
  character(len=*),      intent(in)  :: filename
! maximum length of a line
  integer(kind=i4b),     intent(in)  :: nchars
! Mpi 
  integer(kind=i4b),     intent(in)  :: my_rank 
! number of lines in input file
  integer(kind=i4b),     intent(out) :: nlines
! ----------------------------------------------------------------------
  logical, parameter    :: loutfile = .false.
  integer(kind=i4b)     :: ios, iline
  character(len=nchars) :: line


! check how many lines the input file has
  open(file=filename,unit=iomain,status='old')
  if (loutfile) open(file='out_'//filename,unit=iofile,status='replace')
  nlines = 0
  do
    read(iomain,'(a)',iostat=ios) line
    if (ios /= 0) exit
    nlines = nlines + 1
    if (my_rank == 0) then
      if (loutfile) write(iofile,'(a)') trim(line)
    end if ! my_rank
  end do
  close(iomain)
  if (loutfile) close(iofile)
  if (my_rank == 0) then
    write(*,'(/" Nmber of lines in input=",i8)') nlines
  end if ! my_rank
  if (nlines == 0) stop 'load_input: input file has zero lines!'
! copy input file to memory
  open(file=filename,unit=iomain,status='old')
  if (allocated(inputfile)) deallocate(inputfile)
  allocate(inputfile(nlines))
  do iline=1,nlines
    read(iomain,'(a)',iostat=ios) line
    inputfile(iline) = line
  end do
  close(iomain)
! All done!
  end subroutine load_input
