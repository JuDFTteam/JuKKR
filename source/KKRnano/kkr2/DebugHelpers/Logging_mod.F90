!> A simple logging module.
!> Author: Elias Rabel
!> Use logging macros provided in logging_macros.h instead of
!> the routines provided here. This allows for compilation
!> with or without logging support.
!>
!> Usage:
!> #include 'logging_macros.h'
!> Put macro USE_LOGGING_MOD before 'implicit none'
!> Use OPENLOG once in the program
!> USE WRITELOG to write to log.
!> Use CLOSELOG once when done
!>
!> example
!> WRITELOG(2, *) "Logging ..."
!>          |  |
!>   loglevel  Fortran format string

module Logging_mod
  implicit none
  private
  
  public :: openLogfile, closeLogfile, checkLog
  integer, public, parameter :: LOGFILEHANDLE = 121
  
  integer, private, protected :: logging_level = 0
  logical, private, protected :: log_created = .false.

  CONTAINS

  !----------------------------------------------------------------------------
  !> Creates a logfile numbered by 'rank' and sets 'loglevel'
  !> @param loglevel determines amount of logging 0 (no logging) >0 write log
  subroutine openLogfile(rank, loglevel)
    integer, intent(in) :: rank
    integer, intent(in) :: loglevel

    character(len=6) :: str

    write (str, '(I6.6)') rank

    logging_level = loglevel
    log_created = .false.

    if (loglevel > 0 .and. (log_created .eqv. .false.)) then
      open(LOGFILEHANDLE, file='log.' // str, form='formatted', &
           status='replace')

      log_created = .true.
    end if
  end subroutine

  !----------------------------------------------------------------------------
  !> Closes the logfile
  subroutine closeLogfile()
    if (log_created .eqv. .true.) then
      close(LOGFILEHANDLE)
      log_created = .false.
    end if
  end subroutine

  !----------------------------------------------------------------------------
  !> Checks if a log entry should be written according to given level
  !> Returns true if logging_level >= level
  logical function checkLog(level)
    integer, intent(in) :: level
    checkLog = .false.
    if (logging_level >= level .and. log_created .eqv. .true.) then
      checkLog = .true.
    end if
  end function
  
end module Logging_mod
