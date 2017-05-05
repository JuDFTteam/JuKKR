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
    logical(kind=1)  :: running = .false.
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


  subroutine outTimeStats(self, name, iter, unit)
    type(TimerMpi), intent(inout) :: self !> intent(inout) as this stops the timer
    character(len=*), intent(in)  :: name ! which part of the code was timed
    integer, intent(in), optional :: iter !> iteration number
    integer, intent(in), optional :: unit !> output unit number

    character(len=96) :: stats_string

    integer :: iou
    iou = 2; if(present(unit)) iou = unit
!   if (self%running) call stopTimer(self)
!   if (self%running) write(iou,*) "Warning! Tried to get timer stats when still running ",nowTime - self%start_time ! warning

    stats_string = getTimerStats(self)
    if (present(iter)) then
      write(iou, fmt="('iter:',i4,2x,9a)") iter, name, '  ', trim(stats_string)
    else
      write(iou, fmt="(2x,9a)") name, '  ', trim(stats_string)
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
  subroutine startTimer(self, startTime)
    type(TimerMpi), intent(inout) :: self
    double precision, intent(out), optional :: startTime

    double precision :: nowTime
    nowTime = now()
    if (self%running) then
      write(2,*) "Warning! Tried to start a timer twice at ",nowTime - self%start_time ! warning
    else
      self%running = .true.
      self%start_time = nowTime
    endif
    if (present(startTime)) startTime = nowTime

  endsubroutine ! start or restart

  !----------------------------------------------------------------------------
  subroutine stopTimer(self, stopTime)
    type(TimerMpi), intent(inout) :: self
    double precision, intent(out), optional :: stopTime

    double precision :: nowTime
    nowTime = now()
    if (self%running) then
      self%running = .false.
      call add2TimerStats(self, duration = nowTime - self%start_time)
    else
      write(2,*) "Warning! Tried to stop a timer twice at ",nowTime - self%start_time ! warning
    endif
    if (present(stopTime)) stopTime = nowTime

  endsubroutine ! stop

  !----------------------------------------------------------------------------
  subroutine add2TimerStats(self, duration)
    type(TimerMpi), intent(inout) :: self
    double precision, intent(in) :: duration

    ! add to stats
    self%num_samples = self%num_samples + 1
    self%sum_samples = self%sum_samples + duration
    self%sum_squares = self%sum_squares + duration*duration
    self%min_sample = min(self%min_sample, duration)
    self%max_sample = max(self%max_sample, duration)

  endsubroutine ! add

  !---------------------------------------------------------------------------
  !> Returns elapsed time since last tick or start. Aggregate time interval
  double precision function getTime(self)
    type(TimerMpi), intent(inout) :: self

    double precision :: nowTime

    if (self%running) then
      nowTime = now()
      getTime = nowTime - self%start_time
      call add2TimerStats(self, duration = getTime)
      self%start_time = nowTime ! restart
    else
      getTime = 0.d0
    endif

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
        num,' * ',mean,' +/- ',deviation,' [',self%min_sample,', ',self%max_sample,'] sec'
    case default ; str = '<no timing samples>'
    endselect ! num
  endfunction ! getTimerStats

  !---------------------------------------------------------------------------
  !> This is the only spot where we control the interface towards the MPI library
  double precision function now()
    double precision, external :: MPI_Wtime
    now = MPI_Wtime()
  endfunction ! now

endmodule ! TimerMpi_mod
