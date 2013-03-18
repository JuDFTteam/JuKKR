
!================================================================
!> Collects contributions to the rms errors from all sites and
!> prints rms errors.
!
!> output of a) rms-error
!
! called by main2
!================================================================

subroutine RMSOUT_com(RMSAVQ,RMSAVM,ITER,NSPIN,NAEZ, &
MYLRANK, communicator)

  implicit none

  double precision RMSAVM,RMSAVQ

  integer ITER,NSPIN,NAEZ
  !     ..
  !     .. MPI variables ..
  !     .. L-MPI ..
  integer      MYLRANK, communicator

  double precision RMSQ, RMSM

  call allreduceRMS_com(RMSQ, RMSM, RMSAVQ,RMSAVM,NAEZ, communicator)

  ! ================== MYRANK.EQ.0 =======================================
  if(MYLRANK.eq.0) then

    call printRMSerror(ITER, NSPIN, RMSM, RMSQ)

  end if
  ! ============= MYRANK.EQ.0 ============================================

end

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Communicate the contributions of the rms error of charge and magnetisation
!> density and return results in RMSQ and RMSM.
subroutine allreduceRMS_com(RMSQ, RMSM, &  ! output
RMSAVQ_local,RMSAVM_local,NAEZ, & !in
communicator) !in

  implicit none

  INCLUDE 'mpif.h'

  double precision RMSAVM_local,RMSAVQ_local, &
  WORK1(2),WORK2(2)

  integer NAEZ


  !     .. local scalars ..
  double precision RMSQ,RMSM
  !     ..
  !     .. MPI variables ..
  integer communicator
  integer ierr

  !****************************************************** MPI COLLECT DATA

  WORK1(1) = RMSAVQ_local
  WORK1(2) = RMSAVM_local

  call MPI_ALLREDUCE(WORK1,WORK2,2, &
  MPI_DOUBLE_PRECISION,MPI_SUM,communicator, &
  IERR)

  RMSQ = SQRT(WORK2(1)/NAEZ)
  RMSM = SQRT(WORK2(2)/NAEZ)

  !****************************************************** MPI COLLECT DATA
  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! ======================================================================

end

!------------------------------------------------------------------------------
!> Output rms errors to screen.
subroutine printRMSerror(ITER, NSPIN, RMSM, RMSQ)
  implicit none
  integer :: ITER
  integer :: NSPIN
  double precision :: RMSM
  double precision :: RMSQ

  write(6,'(79(1H-),/)')
  if (NSPIN.eq.2) then
    write (6,fmt=9041) ITER,RMSQ,RMSM
  else
    write (6,fmt=9051) ITER,RMSQ
  end if
  write(6,'(79(1H-))')

9041 format ('      ITERATION',I4,' average rms-error : v+ + v- = ', &
  1p,d11.4,/,39x,' v+ - v- = ',1p,d11.4)
9051 format ('      ITERATION',I4,' average rms-error : v+ + v- = ', &
  1p,d11.4)
end subroutine
