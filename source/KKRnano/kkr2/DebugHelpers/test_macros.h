#ifdef DEBUG1
#define TESTARRAYLOCAL(ARRAY) call testarray(#ARRAY, ARRAY)
#define TESTARRAY(RANK,ARRAY) if (test_getmyrank() == (RANK)) call testarray(#ARRAY // " " // #RANK, ARRAY)
#define USE_ARRAYTEST_MOD use arraytest_mod
#else
#define TESTARRAYLOCAL(ARRAY)
#define TESTARRAY(RANK,ARRAY)
#define USE_ARRAYTEST_MOD
#endif
