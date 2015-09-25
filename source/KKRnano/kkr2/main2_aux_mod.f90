module main2_aux_mod
  implicit none
  private
  public :: printDoubleLineSep, is_abort_by_rank0, writeIterationTimings

  contains

  !----------------------------------------------------------------------------
  !> Print total number of Iterative solver iterations
  subroutine printSolverIterationNumber(ITER, NOITER_ALL)
    integer, intent(in) :: ITER
    integer, intent(in) :: NOITER_ALL

    write(6,'(79(1H=))')
    write(6,'(19X,A,I3,A,I10)') '       ITERATION : ', ITER,' SUM of QMR ',NOITER_ALL
    write(6,'(79(1H=),/)')
  endsubroutine

  !----------------------------------------------------------------------------
  !> Print timing information for SCF iteration.
  subroutine writeIterationTimings(ITER, TIME_I, TIME_S)
    use TimerMpi_mod, only: outtime
    integer, intent(in) :: ITER
    double precision, intent(in) :: TIME_I
    double precision, intent(in) :: TIME_S

    call OUTTIME(.true. ,'end .................', TIME_I, ITER)
    write(6,'(79(1H=))')
    write(2,'(79(1H=))')
    call OUTTIME(.true. ,'finished in .........', TIME_S, ITER)
    write(2,'(79(1H=))')
    write(6,'(79(1H=),/)')
  endsubroutine

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

  endfunction

  !----------------------------------------------------------------------------
  !> Prints a double line separator. (=========)
  !> @param unit_number  optional: write to file 'unit_number' otherwise to
  !>                               stdout
  subroutine printDoubleLineSep(unit_number)
    integer, intent(in), optional :: unit_number

    if (present(unit_number)) then
      write(unit_number,'(79("="))')
    else
      write(*,'(79("="))')
    endif
  endsubroutine

endmodule main2_aux_mod
