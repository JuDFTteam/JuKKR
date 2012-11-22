module muffin_tin_zero_mod
  implicit none

  contains

  !----------------------------------------------------------------------------
  !> Communicate contributions to the muffin-tin-zero (potential)
  !> of each process in 'communicator'
  subroutine allreduceMuffinTinShift_com(communicator, VAV0, VBC, VOL0)
    implicit none

    include 'mpif.h'

    integer :: communicator
    double precision :: VAV0
    double precision :: VBC(2)
    double precision :: VOL0

    !-------------------------------
    integer :: IERR

    double precision :: WORK1(2)
    double precision :: WORK2(2)

    ! Calculate muffin-tin potential shift
    !****************************************************** MPI COLLECT DATA
    WORK1(1) = VAV0
    WORK1(2) = VOL0
    call MPI_ALLREDUCE(WORK1,WORK2,2,MPI_DOUBLE_PRECISION,MPI_SUM, &
    communicator,IERR)
    VAV0 = WORK2(1)
    VOL0 = WORK2(2)
    !****************************************************** MPI COLLECT DATA

    VBC(1) = -VAV0/VOL0
    VBC(2) = VBC(1)
  end subroutine

  !----------------------------------------------------------------------------
  !> Print muffin-tin zero information on screen.
  subroutine printMuffinTinShift(VAV0, VBC, VOL0)
    implicit none
    double precision :: VAV0
    double precision :: VBC(2)
    double precision :: VOL0

    write (6,fmt=9103) VOL0,VAV0,VBC(1)
    write(6,'(79(1H=),/)')

9103 format ('  VOL INT.',F16.9,'  VAV INT.',F16.9,'  VMT ZERO',F16.9)
  end subroutine

end module muffin_tin_zero_mod
