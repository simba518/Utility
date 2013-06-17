#include <iostream>
using namespace std;
#include "MyTestFixture.h"
using namespace TEST;

MyTestFixture::MyTestFixture(){

}

MyTestFixture::~MyTestFixture(){

}

void MyTestFixture::print(const char *test_message){
  if (test_message)
	cout<< endl << "Test: " << test_message << endl;
}
