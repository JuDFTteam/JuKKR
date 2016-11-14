
!================================================================
!> collects contributions to the rms errors from all sites and
!> prints rms errors.
!
!> output of a) rms-error
!
! called by main2
!================================================================

!> rmsavq, rmsavm: on input: local rms erros, on output: global rms
!> errors.
subroutine rmsout_com(rmsavq, rmsavm, iter, nspin, naez, myrank, communicator)
  implicit none

  double precision, intent(inout) :: rmsavm, rmsavq

  integer, intent(in) :: iter, nspin, naez
  integer, intent(in) :: myrank, communicator

  double precision rmsq, rmsm

  call allreducerms_com(rmsq, rmsm, rmsavq, rmsavm, naez, communicator)

  if (myrank == 0) call printrmserror(iter, nspin, rmsm, rmsq)

  rmsavq = rmsq
  rmsavm = rmsm

endsubroutine ! rmsout_com

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> communicate the contributions of the rms error of charge and magnetisation
!> density and return results in rmsq and rmsm.
subroutine allreducerms_com(rmsq, rmsm, rmsavq_local, rmsavm_local, naez, communicator)
  implicit none
  include 'mpif.h'

  double precision, intent(in) :: rmsavm_local, rmsavq_local
  double precision, intent(in) :: rmsq, rmsm
  integer, intent(in) :: naez
  integer, intent(in) :: communicator

  !     .. local scalars ..
  double precision :: send(2), recv(2)
  integer :: ierr

  send(:) = [rmsavq_local, rmsavm_local]

  call MPI_Allreduce(send, recv, 2, mpi_double_precision, mpi_sum, communicator, ierr)

  rmsq = sqrt(recv(1)/naez)
  rmsm = sqrt(recv(2)/naez)

endsubroutine ! allreducerms_com

!------------------------------------------------------------------------------
!> output rms errors to screen.
subroutine printrmserror(iter, nspin, rmsm, rmsq)
  implicit none
  integer, intent(in) :: iter, nspin
  double precision, intent(in) :: rmsm, rmsq

  write(6,'(79(1h-),/)')
  write (6,fmt="('      ITERATION',i4,' average rms-error : v+ + v- = ',1p,d11.4)") iter, rmsq
  if (nspin == 2) &
  write (6,fmt="(39x,' v+ - v- = ',1p,d11.4)") rmsm
  write(6,'(79(1h-))')

endsubroutine ! print
