!> A simple logging module.
!> Author: Elias Rabel
!> Use these logging macros instead of
!> the routines provided in Logging_mod.F90.
!> This allows for compilation with or without logging support.
!>
!> Usage:
!> #include "logging_macros.h"
!> Put macro USE_LOGGING_MOD before 'implicit none'
!> Use OPENLOG once in the program
!> USE WRITELOG to write to log.
!> Use CLOSELOG once when done
!>
!> example
!> WRITELOG(2, *) "Logging ...", xy
!>          |  |
!>   loglevel  Fortran format string
!>
!> If no logging is desired:
!> *) When compiled with macro NOLOGGING defined, then logging code
!> is not included in source.
!> *) Or: when using OPENLOG, set loglevel to 0 - logging code is
!> included - but output is deactivated and no logfile created.

#ifndef LOGGING_MACROS_H_
#define LOGGING_MACROS_H_


#ifndef NOLOGGING
#define USE_LOGGING_MOD use Logging_mod
#define OPENLOG(RANK, LOGLEVEL) call openLogfile( (RANK), (LOGLEVEL) )
#define WRITELOG(LEVEL, FORMAT) if (checkLog((LEVEL))) write(LOGFILEHANDLE, FORMAT)
#define CLOSELOG call closeLogfile()
#else

#define USE_LOGGING_MOD
#define OPENLOG(RANK, LOGLEVEL) continue

! unfortunately a bad trick has to be used
! the effect is to ignore the rest of the line, if it is formatted
! as to be accepted by a write statement
#define WRITELOG(LEVEL, FORMAT) if (.false.) write(*,*)

#define CLOSELOG continue

#endif

#endif
