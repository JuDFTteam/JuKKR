!> TESTARRAYLOG(LOGLEVEL, ARRAY) Writes some diagnostic values calculated from ARRAY to
!> logfile if the loglevel set is >= LOGLEVEL
!> uses the logging system defined in logging_macros.h and Logging_mod
!> Requirements:
!> at beginning of source file:
!> *) #include 'logging_macros.h'
!> *) #include 'test_array_log.h' (this file)
!> in subroutine where logging should be used
!> *) To use the TESTARRAYLOG macro one has to put USE_ARRAYLOG_MOD before 'implicit none'
!> *) and also the USE_LOGGING_MOD macro
!> *) set up the logging system: see logging_macros.h
!> *) compile with preprocessor macro DEBUG1 set
!> Note: the logging system is implemented as macros to allow for compilation
!> with and without logging code

#ifndef TEST_ARRAY_LOG_H_
#define TEST_ARRAY_LOG_H_

#ifdef NOLOGGING
#undef DEBUG1
#endif

#ifdef DEBUG1
#define TESTARRAYLOG(LOGLEVEL, ARRAY) WRITELOG(LOGLEVEL, *) testarray(__LINE__, #ARRAY, ARRAY), __FILE__
#define USE_ARRAYLOG_MOD use arraytest2_mod
#else
#define TESTARRAYLOG(LOGLEVEL, ARRAY)
#define USE_ARRAYLOG_MOD
#endif

#endif
