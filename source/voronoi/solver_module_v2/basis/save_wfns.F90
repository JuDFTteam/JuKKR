  subroutine save_wfns(ia,il,is)
! Output wavefunction for projection to file
  use global

  implicit none

! --> susc atom
  integer(kind=i4b), intent(in) :: ia
! --> angular momentum
  integer(kind=i4b), intent(in) :: il
! --> spin channel
  integer(kind=i4b), intent(in) :: is
! -----------------------------------------------------------------
  integer(kind=i4b) :: i, nb
  character*18      :: filename

  nb = iwsusc(il,is,ia)
  if (nb == 0) return
  write(filename,'("ia",i4.4,"is",i2.2,"il",i2.2,".wfn")') ia,is,il
  open(file=filename,unit=iofile)
  write(iofile,'("# Susc basis set: ia, iasusc, il, is, nb, nr; then rmesh, dr, dr, phiref")')
  write(iofile,'("# ",6i8)') ia, iasusc(ia), il, is, nb, nrpts(ia)
  do i=1,nrpts(ia)
    write(iofile,'(1000es16.8)') rmesh(i,ia), drmesh(i,ia), drproj(i,ia), phiref(i,1:nb,il,is,ia)
  end do
  close(iofile)
! All done!
  end subroutine save_wfns
