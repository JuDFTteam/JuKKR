#ifdef DEBUG1
#define TESTARRAYLOCAL(ARRAY) call testarray(test_getmyrank(), #ARRAY, ARRAY)
#define TESTARRAY(RANK,ARRAY) if (test_getmyrank() == (RANK)) call testarray((RANK), #ARRAY, ARRAY)
#define USE_ARRAYTEST_MOD use arraytest_mod
#else
#define TESTARRAYLOCAL(ARRAY)
#define TESTARRAY(RANK,ARRAY)
#define USE_ARRAYTEST_MOD
#endif
