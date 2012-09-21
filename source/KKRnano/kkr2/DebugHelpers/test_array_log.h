#ifdef DEBUG1
#define TESTARRAYLOG(LOGLEVEL, ARRAY) WRITELOG(LOGLEVEL, *) testarray(__LINE__, #ARRAY, ARRAY), __FILE__
#define USE_ARRAYLOG_MOD use arraytest2_mod
#else
#define TESTARRAYLOG(LOGLEVEL, ARRAY)
#define USE_ARRAYLOG_MOD
#endif
