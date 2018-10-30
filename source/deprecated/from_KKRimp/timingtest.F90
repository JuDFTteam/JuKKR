
include 'timing.F90'

program timingtest
  use mod_timing
  integer :: test
  real*8 :: time
  call timing_init(0)
  call timing_start('bla')
!   read(*,*) test
  call timing_stop('bla',time)
  write(*,*) 'bla',time


  call timing_start('bla2')
  read(*,*) test
  call timing_stop('bla2',time)
  write(*,*) 'bla2',time



end program timingtest


