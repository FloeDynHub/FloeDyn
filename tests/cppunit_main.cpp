// CppUnit version
#include <cppunit/ui/text/TestRunner.h>
#include "examples/qt/ExampleTestCases.h"
#include <iostream>

int main( int argc, char **argv)
{
    using namespace std;
    cout << "ok" << endl;
    
    CppUnit::TextUi::TestRunner runner;

//  runner.addTest( ExampleTestCases::suite() );
//  runner.run();
  return 0;
}