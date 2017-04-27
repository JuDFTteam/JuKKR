  subroutine close_outsusc()
! This just closes the outsusc.dat file
  use global

  implicit none

  close(iomain)
! All done!
  end subroutine close_outsusc

