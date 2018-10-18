!--------------------------------------------------------------------------------
! Copyright (c) 2018 Forschungszentrum Juelich GmbH, Juelich, Germany
! This file is part of KKRnano and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module main2_aux_mod
  implicit none
  private
  public :: is_abort_by_rank0

  contains

  !-------------------------------------------------------------------------
  !> Checks for abort flag on rank 0.
  !>
  !> If rank 0 passes flag>=1 this function returns .true. on all ranks,
  !> otherwise .false.
  !> The value of flag passed by ranks other than rank 0 is ignored.
  logical function is_abort_by_rank0(flag, communicator)
    integer, intent(in) :: flag
    integer, intent(in) :: communicator

    include 'mpif.h'
    integer :: ierr
    integer :: master_rank
    integer :: stop_integer

    master_rank = 0
    stop_integer = flag

    call MPI_BCAST(stop_integer,1,MPI_INTEGER, master_rank,communicator,ierr)

    is_abort_by_rank0 = (stop_integer >= 1)

  endfunction ! is

endmodule ! main2_aux_mod
