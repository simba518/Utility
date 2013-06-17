#ifndef MYTESTFIXTURE_H
#define MYTESTFIXTURE_H

#include <cppunit/TestFixture.h>
#include <cppunit/Asserter.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCaller.h>
#include <cppunit/extensions/HelperMacros.h>

#include <string>
#include <vector>
#include <assert.h>
#include <fstream>

namespace TEST{

  class MyTestFixture:public CppUnit::TestFixture{

  public:
	MyTestFixture();
	~MyTestFixture();

  protected:
	//recording the test information;
	void print(const char *test_message);

  private:
	std::string warped_path;
  };

}


#endif
