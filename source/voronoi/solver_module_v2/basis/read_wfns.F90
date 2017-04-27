  subroutine read_wfns(ia,il,is)
! Read basis functions from file
  use global

  implicit none

! --> susc atom
  integer(kind=i4b), intent(in) :: ia
! --> angular momentum
  integer(kind=i4b), intent(in) :: il
! --> spin channel
  integer(kind=i4b), intent(in) :: is
! -----------------------------------------------------------------
  integer(kind=i4b) :: i, ia1, ia2, il1, is1, nb
  character*1       :: dummy
  character*18      :: filename
  logical           :: exists

  write(filename,'("ia",i4.4,"is",i2.2,"il",i2.2,".wfn")') ia,is,il
! Where is my mind?
  inquire(file=filename,exist=exists)
  if (.not.exists) then
    write(*,*) "in_wfns: file ",filename," not found!"
    stop
  end if
! Was the potential file read first?
  if (normesh(ia)) stop 'in_wfns: read pot first!'
  write(iodb,*) "Reading ", filename
! Now read the basis set
  open(file=filename,unit=iofile,status='old')
  read(iofile,*) ! header
  read(iofile,*) dummy, ia1, ia2, il1, is1, nb
! The number of basis functions is updated here
  iwsusc(il,is,ia) = nb
  do i=1,nrpts(ia)
    read(iofile,*) rmesh(i,ia), drmesh(i,ia), drproj(i,ia), phiref(i,1:nb,il,is,ia)
  end do
  close(iofile)
! wfn was read in
  nowfns(il,is,ia) = .false.
! All done!
  end subroutine read_wfns
