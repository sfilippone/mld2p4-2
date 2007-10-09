#include <sys/resource.h>
#include <unistd.h>
#include <signal.h> 

#ifdef Add_
#define enablecore  enablecore_
#endif
#ifdef AddDouble_
#define enablecore  enablecore_
#endif
#ifdef UpCase
#define enablecore  ENABLECORE_
#endif


void enablecore()
{
  struct rlimit rlim;
  getrlimit(RLIMIT_CORE, &rlim);
  rlim.rlim_cur = RLIM_INFINITY;
  setrlimit(RLIMIT_CORE, &rlim);
  signal(SIGSEGV, SIG_DFL);
}
