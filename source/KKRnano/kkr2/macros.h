#define warn(unit, message)  call launch_warning(message, unit, __FILE__, __LINE__)
#define die_here(message) call die(message, file=__FILE__, line=__LINE__)
#define here trim(__FILE__-":"-__LINE__)
! #define assert(condition) if((condition)) then ; else ; die_here("assert("#condition") failed!") ; endif
#define assert(condition) if((condition)) then ; else ; die_here("assert(condition) failed!") ; endif

! #define CHECKASSERT(X) if (.not. (X)) then; write(*,*) "ERROR: Check " // #X // " failed. ", __FILE__, __LINE__; STOP; endif
#define CHECKASSERT(X) if (.not. (X)) then; write(*,*) "ERROR: Check failed. ", __FILE__, __LINE__; STOP; endif
