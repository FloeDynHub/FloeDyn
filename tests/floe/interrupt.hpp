#include <signal.h>
#include <atomic>

std::atomic<bool> QUIT{false};    // signal flag

namespace interruption 
{
void got_signal(int)
{
    QUIT = true;
}
}