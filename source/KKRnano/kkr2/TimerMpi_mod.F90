!> @author Elias Rabel, 2012

module TimerMpi_mod
  implicit none
  private
  public :: TimerMpi, stopTimer, startTimer, getTime, createTimer, outTime, outTimeStats
  
  type TimerMpi
    private
    double precision :: start_time = 0.d0
    double precision :: sum_samples = 0.d0
    double precision :: sum_squares = 0.d0
    double precision :: min_sample = 0.d0
    double precision :: max_sample = 0.d0
    integer(kind=4)  :: num_samples = 0
    logical :: running = .false.
  endtype ! TimerMpi

  contains

  subroutine outTime(output, name, time, iter)
!>     every time this subroutine is called it produces one line of
!>     output on the file with the unit-number 2, which consists
!>     of up to 59 characters (from the variable name) and the time.
!>     here the time must be in seconds. in the output the time is
!>     given in seconds, but also in hours minutes and seconds.
!>                                              p.kurz   8.2.96
    logical         , intent(in) :: output ! only true for the master rank
    character(len=*), intent(in) :: name ! which part of the code was timed
    double precision, intent(in) :: time ! in seconds
    integer         , intent(in) :: iter !> iteration number

!   double precision :: seconds
!   integer :: ihours, iminutes, idays
!   character(len=*), parameter :: F8000 = "('iter: ',i3,2x,a,t42,f9.2,' sec = ',i3,' h ',i2,' min ',f5.2,' sec')"
    character(len=*), parameter :: F80 = "('iter:',i4,2x,a,t40,f11.2,' sec')" ! print only in seconds

    if (.not. output) return ! only one processor proceeds ...
!     seconds = time
!     
!     idays = floor(seconds/86400.)
!     seconds = seconds - ihours*3600.
!     
!     ihours = floor(seconds/3600.)
!     seconds = seconds - ihours*3600.
! 
!     iminutes = floor(seconds/60.)
!     seconds = seconds - iminutes*60.

!   output of the results to unit 2 and unit 6
!   write(2, fmt=F8000) iter, name, time, ihours, iminutes, seconds ! write to "time-info" unit
    write(2, fmt=F80)   iter, name, time ! write the seconds to "time-info" unit
!   write(6, fmt=F80)   iter, name, time ! skip hours,minutes,seconds display here
  endsubroutine ! outTime

  
  subroutine outTimeStats(self, name, iter)
    type(TimerMpi), intent(inout) :: self !> intent(inout) as this stops the timer
    character(len=*), intent(in)  :: name ! which part of the code was timed
    integer, intent(in), optional :: iter !> iteration number

    character(len=96) :: stats_string
    
    if (self%running) call stopTimer(self)
    stats_string = getTimerStats(self)
    if (present(iter)) then
      write(2, fmt="('iter:',i4,2x,9a)") iter, name, '  ',stats_string
    else
      write(2, fmt="(2x,9a)") name, '  ', stats_string
    endif
  endsubroutine ! outTime
  
  !----------------------------------------------------------------------------
  subroutine createTimer(self)
    type(TimerMpi), intent(inout) :: self

    self%running = .false. ! new timers are stopped when created
    self%start_time = 0.d0
    
    ! init stats    
    self%num_samples = 0
    self%sum_samples = 0.d0
    self%sum_squares = 0.d0
    self%min_sample =  9.d9
    self%max_sample = -9.d9
  endsubroutine ! create
  
  !----------------------------------------------------------------------------
  subroutine startTimer(self)
    type(TimerMpi), intent(inout) :: self
    
    if (self%running) then
      write(2,*) "Warning! Tried to start a timer twice at ",now() - self%start_time ! warning
    else
      self%start_time = now()
      self%running = .true.
    endif
    
  endsubroutine ! start or restart

  !----------------------------------------------------------------------------
  subroutine stopTimer(self)
    type(TimerMpi), intent(inout) :: self

    double precision :: duration
    
    if (self%running) then
      duration = now() - self%start_time
      self%running = .false.
      ! add to stats
      self%num_samples = self%num_samples + 1
      self%sum_samples = self%sum_samples + duration
      self%sum_squares = self%sum_squares + duration*duration
      self%min_sample = min(self%min_sample, duration)
      self%max_sample = max(self%max_sample, duration)
    else
      write(2,*) "Warning! Tried to stop a timer twice at ",now() - self%start_time ! warning
      return ! do nothing since you cannot stop a timer twice
    endif
    
  endsubroutine ! stop

  !---------------------------------------------------------------------------
  !> Returns elapsed time. Look at stop-watch and return the current aggregate time
  double precision function getTime(self)
    type(TimerMpi), intent(in) :: self
    
    double precision :: duration
    duration = 0.d0
    if (self%running) duration = now() - self%start_time
    getTime = self%sum_samples + duration
  endfunction ! get

  character(len=96) function getTimerStats(self) result(str)
    type(TimerMpi), intent(in) :: self

    double precision :: mean, variance, deviation
    integer :: ios, num
    
    num = self%num_samples
    selectcase(num)
    case (1) ; write(unit=str, fmt='(9(f0.3,a))', iostat=ios) self%sum_samples,' sec'
    case (2:)
      mean = self%sum_samples/dble(num)
      variance = max(0.d0, self%sum_squares/dble(num) - mean*mean)
      deviation = sqrt(variance)
      write(unit=str, fmt='(f0.2,a,i0,9(a,f0.3))', iostat=ios) self%sum_samples,' sec = ',&
        num,' * ',mean,' +/- ',deviation,' min ',self%min_sample,' max ',self%max_sample,' sec'
    case default ; str = '<no timing samples>'
    endselect ! num
  endfunction ! getTimerStats

  !---------------------------------------------------------------------------
  !> This is the only spot where we control the interface towards the MPi library
  double precision function now()
    double precision, external :: MPI_Wtime
    now = MPI_Wtime()
  endfunction ! now

endmodule ! TimerMpi_mod
