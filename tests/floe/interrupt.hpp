#include <signal.h>
#include <atomic>

std::atomic<bool> QUIT{false};    // signal flag

namespace interruption 
{
void got_signal(int signum)
{
    printf("Caught interruption signal %d\n",signum);
    QUIT = true;
}
}