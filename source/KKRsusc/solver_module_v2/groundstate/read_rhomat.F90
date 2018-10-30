  subroutine read_rhomat()
! Read density matrix (for LDA+U)
! Call this before ldau_correction
  use global

  implicit none

  integer(kind=i4b) :: ia, ilms, jlms, ia2, nlms2, i2(2), il, is
  character*14      :: filename
  character*1       :: hash
  logical           :: exists
  real(kind=r8b)    :: x(1000)

! **************
  do ia=1,nasusc
! **************
    write(filename,'("rhomat",i4.4,".dat")') ia
    inquire(file=filename,exist=exists)
    if (.not.exists) cycle
    open(file=filename,unit=iofile,status='old')
    read(iofile,*)
    read(iofile,*) hash, ia2, nlms2, vshift(0:nlmax,1:nsmax,ia)
    do ilms=1,nlms
      read(iofile,*) x(1:2*nlms)
      do jlms=1,nlms
        rhomat(ilms,jlms,ia) = cmplx(x(2*jlms-1),x(2*jlms))
      end do
    end do
!   Ueff from average occupation per l and s
    vshift(:,:,ia) = 0.d0
    do ilms=1,nlms
      i2 = i2lms(:,ilms)
      il = i2lm(2,i2(1)); is = i2(2)
      vshift(il,is,ia) = vshift(il,is,ia) + (ueff(il,ia)-jeff(il,ia))*rhomat(ilms,ilms,ia)/(2*il+1.d0)
    end do
    close(iofile)
! ******
  end do
! ******
! All done!
  end subroutine read_rhomat

