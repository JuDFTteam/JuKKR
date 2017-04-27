  subroutine read_rmesh(ia)
! Read radial mesh and potentials
  use global

  implicit none

! Susceptibility atom
  integer(kind=i4b), intent(in) :: ia
! -----------------------------------------------------------------
  integer(kind=i4b) :: nr, ir, nlmax1
  real(kind=r8b)    :: z
  complex(kind=c8b) :: work(nrmax), norm
  character*1       :: dummy
  character*12      :: filename
  logical           :: exists

  write(filename,'("susc",i4.4,".pot")') ia
! Where is my mind?
  inquire(file=filename,exist=exists)
  if (.not.exists) then
    write(*,*) "read_rmesh: file ",filename," not found!"
    stop
  end if
  write(iodb,*) "Reading ", filename
! Read potential file for susc
  open(file=filename,unit=iofile,status='old')
  read(iofile,*) ! header
  read(iofile,*) dummy, z, nr, nlmax1  ! dummy is the symbol #
! Not sure how useful this is, as rs is not used anywhere so far
  if (nlmax1 > nlmax) then
    write(*,*) "read_rmesh: for ia=",ia," nlmax1 > nlmax"
    stop
  end if
! Direction of magnetization
  read(iofile,*) dummy, magdir(:,ia)
! Number of radial points and atomic number
  nrpts(ia) = nr
  zat(ia) = z
  write(iodb,'("nr,z,magdir=",i8,f6.1,3f10.6)') nrpts(ia), zat(ia), magdir(:,ia)
! r, r^s, dr for integration, scalar and magnetic potentials
  do ir=1,nr
    read(iofile,*) rmesh(ir,ia), rsmesh(ir,0:nlmax,ia), drmesh(ir,ia), vr(ir,ia), br(ir,ia), nrc(ir,ia), mrc(ir,ia)
!    write(iodb,'(100es16.8)') rmesh(ir,ia), rsmesh(ir,0:nlmax,ia), drmesh(ir,ia), vr(ir,ia), br(ir,ia), nrc(ir,ia), mrc(ir,ia)
  end do
  close(iofile)
! Mesh and potentials were read
  normesh(ia) = .false.
! All done!
  end subroutine read_rmesh
