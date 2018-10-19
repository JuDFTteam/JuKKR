!#define test
!------------------------------------------------------------------------------------
!> Summary: Container for several routines regarding checksum verification 
!> Author: Philipp Ruessmann
!> The idea of the routine seems to be to verify the checksum of potential files
!> and shapefunction files for unientional corruption. It also contains helper
!> routines for the broadcast
!------------------------------------------------------------------------------------
!> @note  Uncomment the line `#define test` and compile with -D for testing
!> @endnote
!------------------------------------------------------------------------------------
module mod_md5sums

  implicit none

  private
#ifdef CPP_MPI
  public :: md5sum_potential, md5sum_shapefun, lmd5sum_pot, lmd5sum_shape, get_md5sums, mympi_bcast_md5sums
#else
  public :: md5sum_potential, md5sum_shapefun, lmd5sum_pot, lmd5sum_shape, get_md5sums
#endif

  integer :: lmd5sum_pot, lmd5sum_shape
  character (len=:), allocatable :: md5sum_potential, md5sum_shapefun

contains

  !-------------------------------------------------------------------------------  
  !> Summary: Perform the md5 checksum for the potential and shapefunction files
  !> Author: Philipp Ruessmann 
  !> Category: potential, shape-functions, KKRhost
  !> Deprecated: False 
  !> Generates the checksum of the potential and shapefunction files as obtained
  !> via the `md5sum` bash command.
  !-------------------------------------------------------------------------------  
  subroutine get_md5sums(ins, filename_pot, filename_shape)
    implicit none
    integer, intent (in) :: ins !! 0 (MT), 1(ASA), 2(Full Potential)
    character (len=40), intent (in) :: filename_pot   !! Name of the potential file
    character (len=40), intent (in) :: filename_shape !! Name of the shapefunction file
    character (len=100) :: tmp_char
    integer :: istat

    logical :: skipread
#ifndef CPP_NOMD5
    integer :: system
    external :: system
#endif

    skipread = .false.
    ! create md5 checksum for potential
#ifndef CPP_NOMD5
    istat = system('md5 '//trim(filename_pot)//' > tmp_md5sum_pot')
#else
    istat = -1
#endif
    if (istat/=0) then
#ifndef CPP_NOMD5
      istat = system('md5sum '//trim(filename_pot)//' > tmp_md5sum_pot')
#else
      istat = -1
#endif
      if (istat/=0) then
        skipread = .true.
      end if
    end if

    if (skipread) then
      tmp_char = '00000000000000000000000000000000  ' // trim(filename_pot)
    else
      ! read in and then delete temporary file tmp_md5sum_pot
      open (9999, file='tmp_md5sum_pot', status='old', iostat=istat)
      if (istat/=0) stop '[get_md5sums] Error: could not create md5sum for potential file'
      read (9999, '(A100)') tmp_char
      close (9999, status='delete')
    end if                         ! skipread
    lmd5sum_pot = len_trim(tmp_char)
    allocate (character(len=lmd5sum_pot) :: md5sum_potential, stat=istat)
    if (istat/=0) stop '[get_md5sums] Error allocating md5sum_potential'
    md5sum_potential = tmp_char(1:lmd5sum_pot)


    if (ins>0) then
      skipread = .false.
      ! create md5 checksum for shapefun
#ifndef CPP_NOMD5
      istat = system('md5 '//trim(filename_shape)//' > tmp_md5sum_shape')
#else
      istat = -1
#endif
      if (istat/=0) then
#ifndef CPP_NOMD5
        istat = system('md5sum '//trim(filename_shape)//' > tmp_md5sum_shape')
#else
        istat = -1
#endif
        if (istat/=0) then
          skipread = .true.
        end if
      end if

      if (skipread) then
        tmp_char = '00000000000000000000000000000000  ' // trim(filename_shape)
      else
        open (9999, file='tmp_md5sum_shape', status='old', iostat=istat)
        if (istat/=0) stop '[get_md5sums] Error: could not create md5sum for shapefun file'
        read (9999, '(A100)') tmp_char
        close (9999, status='delete')
      end if                       ! skipread
      lmd5sum_shape = len_trim(tmp_char)
      allocate (character(len=lmd5sum_pot) :: md5sum_shapefun, stat=istat)
      if (istat/=0) stop '[get_md5sums] Error allocating md5sum_shapefun'
      md5sum_shapefun = tmp_char(1:lmd5sum_shape)
    end if                         ! ins>0

  end subroutine get_md5sums

#ifdef CPP_MPI
  !-------------------------------------------------------------------------------  
  !> Summary: Broadcast the md5sum checksums to the different nodes 
  !> Author: Philipp Ruessmann 
  !> Category: potential, shape-functions, communication, KKRhost
  !> Deprecated: False 
  !> Broadcast the md5sum checksums to the different nodes
  !-------------------------------------------------------------------------------  
  subroutine mympi_bcast_md5sums(ins, myrank, master)

    use :: mpi
    implicit none
    integer, intent (in) :: ins     !! 0 (MT), 1(ASA), 2(Full Potential)
    integer, intent (in) :: master  !! Index of master process
    integer, intent (in) :: myrank  !! Index of current process
    integer :: ierr

    ! broadcast and allocate checksum for potential
    call mpi_bcast(lmd5sum_pot, 1, mpi_integer, master, mpi_comm_world, ierr)
    if (ierr/=mpi_success) stop '[myMPI_Bcast_md5sums] error broadcasting lmd5sum_pot in main_all'

    if (myrank/=master) then
      allocate (character(len=lmd5sum_pot) :: md5sum_potential, stat=ierr)
      if (ierr/=0) stop '[myMPI_Bcast_md5sums] error allocating md5sum_potential'
    end if

    ! broadcast checksum
    call mpi_bcast(md5sum_potential, lmd5sum_pot, mpi_character, master, mpi_comm_world, ierr)
    if (ierr/=mpi_success) stop '[myMPI_Bcast_md5sums] error broadcasting md5sum_potential in main_all'

    ! do the same for shapefunction
    if (ins>0) then

      call mpi_bcast(lmd5sum_shape, 1, mpi_integer, master, mpi_comm_world, ierr)
      if (ierr/=mpi_success) stop '[myMPI_Bcast_md5sums] error broadcasting lmd5sum_shape in main_all'

      if (myrank/=master) then
        allocate (character(len=lmd5sum_shape) :: md5sum_shapefun, stat=ierr)
        if (ierr/=0) stop '[myMPI_Bcast_md5sums] error allocating md5sum_potential'
      end if
      call mpi_bcast(md5sum_shapefun, lmd5sum_shape, mpi_character, master, mpi_comm_world, ierr)
      if (ierr/=mpi_success) stop '[myMPI_Bcast_md5sums] error broadcasting md5sum_shapefun in main_all'

    end if


  end subroutine mympi_bcast_md5sums
#endif

end module mod_md5sums



! TESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTEST
#ifdef test
!-------------------------------------------------------------------------------  
!> Summary: Test the md5sums routine 
!> Author: Philipp Ruessmann 
!> Category: potential, shape-functions, unit-test, KKRhost
!> Deprecated: False 
!> Test the md5sums routine
!-------------------------------------------------------------------------------  
program test

  use :: mod_md5sums
  implicit none

  call get_md5sums(1, 'potential', 'shapefun')

  write (*, *) 'md5sum for potential file:'
  write (*, '(A)') md5sum_potential

  write (*, *) 'md5sum for shapefun file:'
  write (*, '(A)') md5sum_shapefun

end program test
#endif
! TESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTEST
