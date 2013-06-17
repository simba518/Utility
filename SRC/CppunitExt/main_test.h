#include <cppunit/TestResult.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

#include <iostream>
using namespace std;

namespace TEST{

  /** 
   * the actual main function of the test application.
   * it use argv to control which test suits to be runned.
   * should be called in the main function of the test application.
   * 
   * @param argc the length of argv, if argc==1 ,all test suites will be runned.
   * @param argv argv[0] contain the application path, and others(argv[n],n >= 1)
   * including the test suites' names to be test.
   * @return 0
   */
  int main_test(int argc, char *argv[]){

	CppUnit::TestFactoryRegistry &registry = CppUnit::TestFactoryRegistry::getRegistry("No Test");

	if(argc <= 1){
	  registry.registerFactory( &CppUnit::TestFactoryRegistry::getRegistry("All Tests") );
	}else{
	  for(int i = 1; i < argc; i ++){
		registry.registerFactory( &CppUnit::TestFactoryRegistry::getRegistry(argv[i]) );
	  }
	}

	CppUnit::TextUi::TestRunner runner;
	runner.addTest(registry.makeTest());
	runner.run();
  
	return 0;
  }

}
