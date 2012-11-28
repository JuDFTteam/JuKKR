#ifdef DEBUG1
#define TESTARRAYLOCAL(ARRAY) call testarray(test_getmyrank(), #ARRAY, ARRAY)
#define TESTARRAY(RANK,ARRAY) if (test_getmyrank() == (RANK)) call testarray((RANK), #ARRAY, ARRAY)
#define USE_ARRAYTEST_MOD use arraytest_mod
#else
#define TESTARRAYLOCAL(ARRAY)
#define TESTARRAY(RANK,ARRAY)
#define USE_ARRAYTEST_MOD
#endif

#ifndef NDEBUG
#define ASSERT(CONDITION) if (.not. (CONDITION)) then; write(*,*) "Assertion ", #CONDITION, " failed: ",  __FILE__, __LINE__; endif
#else
#define ASSERT(CONDITION)
#endif

#define CHECKASSERT(X) if (.not. (X)) then; write(*,*) "ERROR: Check " // #X // " failed. ", __FILE__, __LINE__; endif
