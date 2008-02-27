#include <sys/resource.h>
#include <unistd.h>
#include <signal.h> 

#ifdef  LowerUnderscore
#define enablecore  enablecore_
#endif
#ifdef  LowerDoubleUnderscore
#define enablecore  enablecore_
#endif
#ifdef  LowerCase
#define enablecore  enablecore
#endif
#ifdef  UpperUnderscore
#define enablecore  ENABLECORE_
#endif
#ifdef  UpperDoubleUnderscore
#define enablecore  ENABLECORE_
#endif
#ifdef  UpperCase
#define enablecore  ENABLECORE
#endif


void enablecore()
{
  struct rlimit rlim;
  getrlimit(RLIMIT_CORE, &rlim);
  rlim.rlim_cur = RLIM_INFINITY;
  setrlimit(RLIMIT_CORE, &rlim);
  signal(SIGSEGV, SIG_DFL);
}
