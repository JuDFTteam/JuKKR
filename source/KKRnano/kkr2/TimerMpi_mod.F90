!> @author Elias Rabel, 2012

module TimerMpi_mod
  implicit none
  private
  public :: TimerMpi, resetTimer, stopTimer, resumeTimer, getElapsedTime, getInitialTime
  public :: outTime
  
  type TimerMpi
    private
    double precision :: initial_time = 0.d0
    double precision :: intermediate_time = 0.d0
    double precision :: elapsed_time = 0.d0
    logical :: stopped = .true.
  endtype ! TimerMpi

  contains

  subroutine outTime(output, name, time_i, iter)
!>     every time this subroutine is called it produces one line of
!>     output on the file with the unit-number 2, which consists
!>     of up to 59 characters (from the variable name) and the time.
!>     here the time must be in seconds. in the output the time is
!>     given in seconds, but also in hours minutes and seconds.
!>                                              p.kurz   8.2.96
    logical         , intent(in) :: output
    character(len=*), intent(in) :: name
    double precision, intent(in) :: time_i
    integer         , intent(in) :: iter

    double precision :: rest, seconds
    integer :: ihours, iminutes
    character(len=*), parameter :: F8000 = "('iter: ',i3,2x,a,t42,f9.2,' sec = ',i3,' h ',i2,' min ',f5.2,' sec')"
    character(len=*), parameter :: F80 = "('iter: ',i3,2x,a,t42,f9.2,' sec')" ! print only in seconds

    if (output) then ! only one processor proceeds ...
      ! calculate time in hours, minutes and seconds
      rest = time_i
      ihours = int(rest/3600.0)
      rest = rest - real(ihours)*3600

      iminutes = int(rest/60.0)
      seconds = rest - real(iminutes)*60

      ! output of the results to unit 2 and unit 6
!     write(2, fmt=F8000) iter,name,time_i,ihours,iminutes,seconds ! write to "time-info" unit
      write(2, fmt=F80)   iter,name,time_i ! write the seconds to "time-info" unit
!     write(6, fmt=F80)   iter,name,time_i ! skip hours,minutes,seconds display here
    endif ! output
  endsubroutine ! outTime
  
  
  !----------------------------------------------------------------------------
  subroutine resetTimer(timer)
    type(TimerMpi), intent(inout) :: timer

    timer%initial_time = wall_clock_time()
    timer%intermediate_time = timer%initial_time
    timer%elapsed_time = 0.d0
    timer%stopped = .false.
  endsubroutine ! reset

  !----------------------------------------------------------------------------
  subroutine stopTimer(timer)
    type(TimerMpi), intent(inout) :: timer

    double precision :: duration
    if (.not. timer%stopped) then
      duration = wall_clock_time() - timer%intermediate_time
      timer%elapsed_time = timer%elapsed_time + duration
      timer%stopped = .true.
    endif
  endsubroutine ! stop

  !----------------------------------------------------------------------------
  subroutine resumeTimer(timer)
    type(TimerMpi), intent(inout) :: timer

    if (timer%stopped) then
      timer%intermediate_time = wall_clock_time()
      timer%stopped = .false.
    endif
  endsubroutine ! resume

  !---------------------------------------------------------------------------
  !> Returns elapsed time.
  !> Look at stop-watch and return time
  double precision function getElapsedTime(timer)
    type(TimerMpi), intent(in) :: timer
    
    double precision :: duration
    duration = 0.d0
    if (.not. timer%stopped) duration = wall_clock_time() - timer%intermediate_time
    getElapsedTime = timer%elapsed_time + duration
  endfunction ! get

  !---------------------------------------------------------------------------
  double precision function getInitialTime(timer)
    type(TimerMpi), intent(in) :: timer
    
    getInitialTime = timer%initial_time
  endfunction ! get
  
  !---------------------------------------------------------------------------
  !> This is the only spot where we control the interface towards the MPi library
  double precision function wall_clock_time() result(now)
    include 'mpif.h'
    now = MPI_Wtime()
  endfunction ! wtime

endmodule ! TimerMpi_mod
