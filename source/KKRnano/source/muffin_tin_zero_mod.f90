!--------------------------------------------------------------------------------
! Copyright (c) 2018 Forschungszentrum Juelich GmbH, Juelich, Germany
! This file is part of KKRnano and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module muffin_tin_zero_mod
  implicit none
  private
  public :: allreduceMuffinTinShift_com, printMuffinTinShift
  
  contains

  !----------------------------------------------------------------------------
  !> Communicate contributions to the muffin-tin-zero (potential)
  !> of each process in 'communicator'
  subroutine allreduceMuffinTinShift_com(communicator, VAV0, VBC, VOL0)
    include 'mpif.h'

    integer, intent(in) :: communicator
    double precision, intent(inout) :: VAV0
    double precision, intent(out)   :: VBC(2)
    double precision, intent(inout) :: VOL0

    !-------------------------------
    integer :: IERR

    double precision :: WORK1(2)
    double precision :: WORK2(2)

    ! Calculate muffin-tin potential shift
    !****************************************************** MPI COLLECT DATA
    WORK1(1) = VAV0
    WORK1(2) = VOL0
    call MPI_ALLREDUCE(WORK1,WORK2,2,MPI_DOUBLE_PRECISION,MPI_SUM, communicator,IERR)
    VAV0 = WORK2(1)
    VOL0 = WORK2(2)
    !****************************************************** MPI COLLECT DATA

    VBC(1) = -VAV0/VOL0
    VBC(2) = VBC(1)
  end subroutine

  !----------------------------------------------------------------------------
  !> Print muffin-tin zero information on screen.
  subroutine printMuffinTinShift(VAV0, VBC, VOL0)
    double precision, intent(in) :: VAV0
    double precision, intent(in) :: VBC(2)
    double precision, intent(in) :: VOL0

    write(6,fmt="('  VOL INT.',F16.9,'  VAV INT.',F16.9,'  VMT ZERO',F16.9)") VOL0,VAV0,VBC(1)
    write(6,'(79(1H=),/)')
  end subroutine

end module muffin_tin_zero_mod
