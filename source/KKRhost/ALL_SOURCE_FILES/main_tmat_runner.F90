!------------------------------------------------------------------------------------
!> Summary: Dummy version of `main_all` which is used to test the `tmat_newsolver`
!> Author: Philipp Ruessmann
!> Category: unit-test 
!> Deprecated: False 
!> This is dummy version of `main_all`, used to execute first `main0` and then 
!> `main1a_dummy` in which `tmat_newsolver` is called
!------------------------------------------------------------------------------------
program tmat_runner

  use :: constants
  use :: mod_profiling
  use :: mod_main0
  use :: mod_types
  use :: mod_timing
  use :: mod_md5sums
  use :: memoryhandling
  use :: mod_version_info
  use :: global_variables
  use :: mod_datatypes, only: dp
  use :: mod_mympi, only: mympi_init, myrank, nranks, master, mpiatom, mpiadapt
  use :: mod_save_wavefun, only: t_wavefunctions

  implicit none

  integer :: ierr, i_stat, i_all
  character (len=3) :: ctemp                               ! name for output file


  call mympi_init()
  write (ctemp, '(I03.3)') myrank
  call construct_serialnr()
  call memocc(0, 0, 'count', 'start')
  call timing_init(myrank)
  call timing_start('main0')

  ! open output files
  open (1337, file='output.'//trim(ctemp)//'.txt')
  call version_print_header(1337)

  t_inc%i_write = 1
  call main0()
  call timing_stop('main0')

  t_inc%i_write = 1
  t_inc%i_time = 1

  call timing_start('main1a')
  call main1a_dummy()
  call timing_stop('main1a')

  ! delete temporary files
  open (69, file='abvmad.unformatted')
  close (69, status='delete')

  ! print memory report to stdout
  if (t_inc%i_write>0) call memocc(0, 0, 'count', 'stop')

end program tmat_runner

