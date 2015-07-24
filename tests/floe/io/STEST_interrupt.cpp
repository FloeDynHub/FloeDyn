#include <iostream>
#include <signal.h>
// #include <cstring>
#include <atomic>
#include <chrono>

std::atomic<bool> quit{false};    // signal flag

void got_signal(int)
{
    quit = true;
}

class Foo
{
public:
    ~Foo() { std::cout << "destructor\n"; }
};

int main(void)
{
    using namespace std;
    struct sigaction sa;
    memset( &sa, 0, sizeof(sa) );
    sa.sa_handler = got_signal;
    sigfillset(&sa.sa_mask);
    sigaction(SIGINT,&sa,NULL);

    Foo foo;    // needs destruction before exit
    int i = 0;
    auto t_start = chrono::high_resolution_clock::now();
    while (i<1e7)
    {
        // do real work here...
        i++;
        if( quit ) break;    // exit normally after SIGINT
    }
    auto t_end = chrono::high_resolution_clock::now();
    cout << "with break : " << chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;

    t_start = chrono::high_resolution_clock::now();
    while (i<1e7)
    {
        // do real work here...
        i++;
    }
    t_end = chrono::high_resolution_clock::now();
    cout << "without break : " << chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;
    return 0;
}