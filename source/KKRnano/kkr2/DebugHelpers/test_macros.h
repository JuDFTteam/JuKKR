
#ifndef NDEBUG
! #define ASSERT(CONDITION) if (.not. (CONDITION)) then; write(*,*) "Assertion ", #CONDITION, " failed: ",  __FILE__, __LINE__; STOP; endif
#define ASSERT(CONDITION) if (.not. (CONDITION)) then; write(*,*) "Assertion failed: ",  __FILE__, __LINE__; STOP; endif
#else
#define ASSERT(CONDITION)
#endif

! #define CHECKASSERT(X) if (.not. (X)) then; write(*,*) "ERROR: Check " // #X // " failed. ", __FILE__, __LINE__; STOP; endif
#define CHECKASSERT(X) if (.not. (X)) then; write(*,*) "ERROR: Check failed. ", __FILE__, __LINE__; STOP; endif
