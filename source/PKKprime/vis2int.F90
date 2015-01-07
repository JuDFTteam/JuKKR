program findvis2int
  use mod_mympi,      only: mympi_init, myrank, nranks, master
  use mod_read,       only: read_kpointsfile_vis, read_kpointsfile_int
  use mod_fermisurf_basic,  only: save_kpointsfile_vis
#ifdef CPP_MPI
  use mpi
#endif

  implicit none


  integer :: ierr, ikp, n_status_out

  !symmetry arrays
  integer :: nsym1, nsym2
  integer, allocatable :: isym1(:), isym2(:)

  !local k-point arrays
  integer :: nkpts_irr, nkpts_int, nkpts_vis
  integer, allocatable :: kpt2irr(:), irr2kpt(:), vis2int(:), vis2int_tmp(:)
  double precision, allocatable :: kpoints_vis(:,:), kpoints_irr(:,:), kpoints_int(:,:), areas(:)

!=================================!
!========== BEGIN CODE ===========!

  !init
#ifdef CPP_MPI
  call MPI_Init( ierr )
#endif
  call mympi_init()

  call read_kpointsfile_vis(nkpts_vis, nkpts_irr, kpoints_irr, nsym1, isym1, kpt2irr, irr2kpt)
  call read_kpointsfile_int(nkpts_int, kpoints_int, areas, nsym2, isym2)

  allocate(kpoints_vis(3,nkpts_vis), STAT=ierr)
  if(ierr/=0) stop 'Problem allocating kpoints_vis'
  kpoints_vis = kpoints_irr(:,kpt2irr)

  if(nsym1/=nsym2)then
   stop 'wrong number of symmetries'
  else
   if(any(isym1/=isym2))then
     stop 'symmetries not compatible'
   end if
  end if

  allocate(vis2int(nkpts_vis), vis2int_tmp(nkpts_irr), STAT=ierr)
  if(ierr/=0) stop 'Problem allocating vis2int'

  write(*,'(A,I0,A)') 'Finding the nearest neighbours to ', nkpts_vis, ' kpoints.'
  n_status_out = nkpts_irr/100

  do ikp=1,nkpts_irr
    call find_closest_kpoint(nkpts_int, kpoints_int, kpoints_irr(:,ikp), vis2int_tmp(ikp))
!   write(*,*) ikp
    if(mod(ikp,n_status_out)==0) write(*,'(I0,A)') nint(float(ikp*100)/float(nkpts_irr)), ' percent done'
  end do!ikp

  if(any(vis2int_tmp==-1)) stop 'something went wrong'

  vis2int = vis2int_tmp(kpt2irr)

  call save_kpointsfile_vis(nkpts_vis, nkpts_irr, kpoints_vis, nsym1, isym1, kpt2irr, irr2kpt, vis2int, 'FSpoints.vis2int.new')

  call MPI_Finalize( ierr )


contains

  subroutine find_closest_kpoint(nkpts, kpoints_arr, refpoint, minindex)
    implicit none
    integer, intent(in) :: nkpts
    double precision, intent(in) :: kpoints_arr(3,nkpts), refpoint(3)
    integer, intent(out) :: minindex

    integer :: ikp2
    double precision :: distvec(nkpts), mindist, dist

    minindex = -1
    mindist = 1d18

    do ikp2=1,nkpts
      dist = sum((kpoints_arr(:,ikp2) - refpoint(:))**2)
      if(dist<mindist)then
        mindist = dist
        minindex = ikp2
      end if!dist<mindist
    end do!ikp

  end subroutine find_closest_kpoint


end program findvis2int
