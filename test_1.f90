program TEST

  use mod_read,       only: read_kpointsfile_vis
  use mod_symmetries, only: pointgrp, rotate_kpoints
  use mod_parutils,   only: parallel_quicksort
  use mod_mympi,      only: mympi_init, myrank, nranks, master
  use mpi
  implicit none

  integer                       :: nkpts1, nkpts
  integer, allocatable          :: sortindex(:)
  double precision, allocatable :: kpoints1(:,:), kpoints(:,:), values(:)

  integer :: nsym
  integer, allocatable :: isym(:)
  double precision :: rotmat(64,3,3)
  character(len=10) :: rotname(64)

  integer :: ierr, ii
  double precision :: bounds(3)

  !initialize MPI
  call MPI_Init ( ierr )
  call mympi_init()

  call pointgrp(rotmat,rotname) 

  call read_kpointsfile_vis(nkpts1, kpoints1, nsym, isym)
  call rotate_kpoints(rotmat, nkpts1, kpoints1, nsym, isym, nkpts, kpoints)



  do ii=1,3
    bounds(ii) = minval(kpoints(ii,:))
  end do

  allocate(values(nkpts), sortindex(nkpts), STAT=ierr)
  if(ierr/=0) stop 'Problem allocating values'
  values = 1d0*(kpoints(1,:)-bounds(1)) + 3d0*(kpoints(2,:)-bounds(1)) + 5d0*(kpoints(3,:)-bounds(1))

  if(myrank==master) write(*,*) 'start'
  call parallel_quicksort(nranks, myrank, master, nkpts, values, sortindex)
  if(myrank==master) write(*,*) 'stop'

  call MPI_Finalize(ierr)

end program TEST
