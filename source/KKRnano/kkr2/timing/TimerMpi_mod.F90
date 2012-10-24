!> @author Elias Rabel, 2012

module TimerMpi_mod
  implicit none

  type TimerMpi
    private
    double precision :: initial_time = 0.0d0
    double precision :: intermediate_time = 0.0d0
    double precision :: elapsed_time = 0.0d0
    logical :: stopped = .true.
  end type

  CONTAINS

  !----------------------------------------------------------------------------
  subroutine resetTimer(timer)
    implicit none
    type (TimerMpi), intent(inout) :: timer
    include 'mpif.h'

    timer%initial_time = MPI_WTIME()
    timer%intermediate_time = timer%initial_time
    timer%elapsed_time = 0.0d0
    timer%stopped = .false.
  end subroutine

  !----------------------------------------------------------------------------
  subroutine stopTimer(timer)
    implicit none
    type (TimerMpi), intent(inout) :: timer
    include 'mpif.h'

    double precision :: the_time

    if (timer%stopped .eqv. .false.) then
      the_time = MPI_WTIME() - timer%intermediate_time
      timer%elapsed_time = timer%elapsed_time + the_time
      timer%stopped = .true.
    end if
  end subroutine

  !----------------------------------------------------------------------------
  subroutine resumeTimer(timer)
    implicit none
    type (TimerMpi), intent(inout) :: timer
    include 'mpif.h'

    if (timer%stopped .eqv. .true.) then
      timer%intermediate_time = MPI_WTIME()
      timer%stopped = .false.
    end if
  end subroutine

  !---------------------------------------------------------------------------
  !> Returns elapsed time.
  !> Look at stop-watch and return time
  double precision function getElapsedTime(timer)
    implicit none
    type (TimerMpi), intent(in) :: timer
    include 'mpif.h'
    double precision :: the_time

    the_time = 0.0d0

    if (timer%stopped .eqv. .false.) then
      the_time = MPI_WTIME() - timer%intermediate_time
    end if

    getElapsedTime = timer%elapsed_time + the_time
  end function

  !---------------------------------------------------------------------------
  double precision function getInitialTime(timer)
    implicit none
    type (TimerMpi), intent(in) :: timer
    getInitialTime = timer%initial_time
  end function

end module TimerMpi_mod
