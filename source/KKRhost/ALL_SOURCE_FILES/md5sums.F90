! uncomment the following line and compile with -D for testing
!#define test


module mod_md5sums

  Use mod_datatypes, Only: dp
  implicit none

  private
#ifdef CPP_MPI
  public md5sum_potential, md5sum_shapefun, lmd5sum_pot, lmd5sum_shape, get_md5sums, myMPI_Bcast_md5sums
#else
  public md5sum_potential, md5sum_shapefun, lmd5sum_pot, lmd5sum_shape, get_md5sums
#endif
 
  integer :: lmd5sum_pot, lmd5sum_shape
  character(len=:), allocatable :: md5sum_potential, md5sum_shapefun
 
  contains
 
  subroutine get_md5sums(ins, filename_pot, filename_shape)
    implicit none
    integer, intent(in) :: ins
    character(len=40), intent(in) :: filename_pot, filename_shape
    character(len=100) :: tmp_char
    integer :: istat

    logical :: reorder, skipread
#ifndef CPP_NOMD5
    integer :: SYSTEM
    external :: SYSTEM
#endif
    
    reorder  = .false.
    skipread = .false.
    ! create md5 checksum for potential
#ifndef CPP_NOMD5
    istat = SYSTEM('md5 '//trim(filename_pot)//' > tmp_md5sum_pot')
#else
    istat=-1
#endif
    if(istat/=0)then
#ifndef CPP_NOMD5
      istat = SYSTEM('md5sum '//trim(filename_pot)//' > tmp_md5sum_pot')
#else
      istat = -1
#endif
      reorder=.true.
      if(istat/=0)then
        skipread=.true.
      end if
    end if

    if(skipread) then
      tmp_char = '00000000000000000000000000000000  '//trim(filename_pot)
    else
      ! read in and then delete temporary file tmp_md5sum_pot
      open(9999, file='tmp_md5sum_pot', status='old', iostat=istat)
      if(istat/=0) stop '[get_md5sums] Error: could not create md5sum for potential file'
      read(9999, '(A100)') tmp_char
      close(9999, status='delete')
    end if!skipread
    lmd5sum_pot = len_trim(tmp_char)
    allocate(character(len=lmd5sum_pot) :: md5sum_potential, stat=istat)
    if(istat/=0) stop '[get_md5sums] Error allocating md5sum_potential'
    md5sum_potential = tmp_char(1:lmd5sum_pot)




    if(ins>0) then
      reorder  = .false.
      skipread = .false.
      ! create md5 checksum for shapefun
#ifndef CPP_NOMD5
      istat = SYSTEM('md5 '//trim(filename_shape)//' > tmp_md5sum_shape')
#else
      istat=-1
#endif
      if(istat/=0)then
#ifndef CPP_NOMD5
        istat = SYSTEM('md5sum '//trim(filename_shape)//' > tmp_md5sum_shape')
#else
        istat=-1
#endif
        reorder=.true.
        if(istat/=0)then
          skipread=.true.
        end if
      end if

      if(skipread) then
        tmp_char = '00000000000000000000000000000000  '//trim(filename_shape)
      else
        open(9999, file='tmp_md5sum_shape', status='old', iostat=istat)
        if(istat/=0) stop '[get_md5sums] Error: could not create md5sum for shapefun file'
        read(9999, '(A100)') tmp_char
        close(9999, status='delete')
      end if!skipread
      lmd5sum_shape = len_trim(tmp_char)
      allocate(character(len=lmd5sum_pot) :: md5sum_shapefun, stat=istat)
      if(istat/=0) stop '[get_md5sums] Error allocating md5sum_shapefun'
      md5sum_shapefun = tmp_char(1:lmd5sum_shape)
    end if!ins>0
  
  end subroutine get_md5sums
  
#ifdef CPP_MPI
  subroutine myMPI_Bcast_md5sums(ins, myrank, master)
    
    use mpi
    implicit none
    integer, intent(in) :: ins, myrank, master
    integer :: ierr
    
    ! broadcast and allocate checksum for potential
    call MPI_Bcast(lmd5sum_pot, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_SUCCESS) stop '[myMPI_Bcast_md5sums] error broadcasting lmd5sum_pot in main_all'
    
    if(myrank/=master) then
      allocate(character(len=lmd5sum_pot) :: md5sum_potential, stat=ierr)
      if(ierr/=0) stop '[myMPI_Bcast_md5sums] error allocating md5sum_potential'
    end if
    
    ! broadcast checksum
    call MPI_Bcast(md5sum_potential, lmd5sum_pot, MPI_CHARACTER, master, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_SUCCESS) stop '[myMPI_Bcast_md5sums] error broadcasting md5sum_potential in main_all'
    
    
    
    ! do the same for shapefunction
    if(ins>0) then
    
      call MPI_Bcast(lmd5sum_shape, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      if(ierr/=MPI_SUCCESS) stop '[myMPI_Bcast_md5sums] error broadcasting lmd5sum_shape in main_all'
    
      if(myrank/=master) then
        allocate(character(len=lmd5sum_shape) :: md5sum_shapefun, stat=ierr)
        if(ierr/=0) stop '[myMPI_Bcast_md5sums] error allocating md5sum_potential'
      end if
      call MPI_Bcast(md5sum_shapefun, lmd5sum_shape, MPI_CHARACTER, master, MPI_COMM_WORLD, ierr)
      if(ierr/=MPI_SUCCESS) stop '[myMPI_Bcast_md5sums] error broadcasting md5sum_shapefun in main_all'
      
    end if
    
    
  end subroutine myMPI_Bcast_md5sums
#endif

end module mod_md5sums



!TESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTEST
#ifdef test
program test

  use mod_md5sums
  implicit none
 
  call get_md5sums(1, 'potential', 'shapefun')
 
  write(*,*) 'md5sum for potential file:'
  write(*,'(A)') md5sum_potential
 
  write(*,*) 'md5sum for shapefun file:'
  write(*,'(A)') md5sum_shapefun

end program test
#endif
!TESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTEST
