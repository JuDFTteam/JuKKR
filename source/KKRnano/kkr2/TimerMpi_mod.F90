!> @author Elias Rabel, 2012

module TimerMpi_mod
  implicit none
  private
  public :: TimerMpi, resetTimer, stopTimer, resumeTimer, getElapsedTime, getInitialTime
  ! new
  public :: outtime
  
  type TimerMpi
    private
    double precision :: initial_time = 0.0d0
    double precision :: intermediate_time = 0.0d0
    double precision :: elapsed_time = 0.0d0
    logical :: stopped = .true.
  end type

  CONTAINS

      subroutine outtime(output, name, time_i, iter)
!>**********************************************************************
!>     every time this subroutine is called it produces one line of
!>     output on the file with the unit-number 2, which consitst
!>     of up to 59 characters (from the variable name) and the time.
!>     here the time must be in seconds. in the output the time is
!>     given in seconds, but also in hours minutes and seconds.
!>                                              p.kurz   8.2.96
!>**********************************************************************
      logical         , intent(in) :: output
      character(len=*), intent(in) :: name
      double precision, intent(in) :: time_i
      integer         , intent(in) :: iter

      double precision :: rest, seconds
      integer :: ihours, iminutes
      character(len=*), parameter :: F8000 = "('iter: ',i3,2x,a,t42,f9.2,' sec = ',i3,' h ',i2,' min ',f5.2,' sec')"

!     only one processor proceeds ...
      if (output) then
!       calculate time in hours, minutes and seconds
        rest = time_i
        ihours = int(rest/3600.0)
        rest = rest - real(ihours)*3600

        iminutes = int(rest/60.0)
        seconds = rest - real(iminutes)*60

!       output of the results to unit 2 and unit 6

        write (2,fmt=F8000) iter,name,time_i,ihours,iminutes,seconds
        write (6,fmt=F8000) iter,name,time_i,ihours,iminutes,seconds
      endif ! output
! 8000 format ('iter: ',i3,2x,a,t42,f9.2,' sec = ',i3,' h ',i2,' min ',f5.2,' sec')
      end subroutine outtime
  
  
  !----------------------------------------------------------------------------
  subroutine resetTimer(timer)
    type (TimerMpi), intent(inout) :: timer

    timer%initial_time = wall_clock_time()
    timer%intermediate_time = timer%initial_time
    timer%elapsed_time = 0.0d0
    timer%stopped = .false.
  end subroutine

  !----------------------------------------------------------------------------
  subroutine stopTimer(timer)
    type (TimerMpi), intent(inout) :: timer

    double precision :: the_time

    if (timer%stopped .eqv. .false.) then
      the_time = wall_clock_time() - timer%intermediate_time
      timer%elapsed_time = timer%elapsed_time + the_time
      timer%stopped = .true.
    end if
  end subroutine

  !----------------------------------------------------------------------------
  subroutine resumeTimer(timer)
    type (TimerMpi), intent(inout) :: timer

    if (timer%stopped .eqv. .true.) then
      timer%intermediate_time = wall_clock_time()
      timer%stopped = .false.
    end if
  end subroutine

  !---------------------------------------------------------------------------
  !> Returns elapsed time.
  !> Look at stop-watch and return time
  double precision function getElapsedTime(timer)
    type (TimerMpi), intent(in) :: timer
    
    double precision :: the_time

    the_time = 0.0d0

    if (timer%stopped .eqv. .false.) then
      the_time = wall_clock_time() - timer%intermediate_time
    end if

    getElapsedTime = timer%elapsed_time + the_time
  end function

  !---------------------------------------------------------------------------
  double precision function getInitialTime(timer)
    type (TimerMpi), intent(in) :: timer
    
    getInitialTime = timer%initial_time
  end function
  
  !---------------------------------------------------------------------------
  !> This is the only spot where we control the interface towards the MPi library
  double precision function wall_clock_time() result(now)
    include 'mpif.h'
    now = MPI_WTIME()
  endfunction

end module TimerMpi_mod
