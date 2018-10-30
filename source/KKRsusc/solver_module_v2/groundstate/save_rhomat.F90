  subroutine save_rhomat()
! Save density matrix (for LDA+U)
! Reduced precision to speed up convergence
  use global

  implicit none

  integer(kind=i4b) :: ia, ilms
  character*14      :: filename

! **************
  do ia=1,nasusc
! **************
!   Write file for rhomat
    write(filename,'("rhomat",i4.4,".dat")') ia
    open(file=filename,unit=iofile,status='replace')
    write(iofile,'("# rhomat: ia, nlms, vshift")')
    write(iofile,'("# ",2i8,100f12.6)') ia, nlms, vshift(0:nlmax,1:nsmax,ia)
    do ilms=1,nlms
      write(iofile,'(100f8.4)') rhomat(ilms,1:nlms,ia)
    end do
    close(iofile)
! ******
  end do
! ******
! All done!
  end subroutine save_rhomat

