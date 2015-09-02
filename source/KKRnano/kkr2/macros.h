#define warn(unit, message)  call launch_warning(message, unit, __FILE__, __LINE__)
#define die_here(message) call die(message, file=__FILE__, line=__LINE__)
#define here trim(__FILE__-":"-__LINE__)
#define assert(condition) if((condition)) then ; else ; die_here("assert("#condition") failed!") ; endif
