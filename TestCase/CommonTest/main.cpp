#define LOG_WARN
#define LOG_ERROR
#define LOG_DEBUG
#define LOG_INFO
#define LOG_TRACE

#include <JsonFilePaser.h>
#include <Timer.h>
#include <IO.h>
#include <assertext.h>
#include <iostream>
#include <Log.h>
#include <MatrixIO.h>
#include <eigen3/Eigen/Dense>
using namespace UTILITY;
using namespace std;
using namespace Eigen;

int main(int argc, char *argv[]){

  TRACE_FUN();
  ERROR_LOG("asd;asjd;fj");
  ERROR_LOG_COND("asd;asjd;fj",false);
  WARN_LOG("asd;asjd;fj");
  DEBUG_LOG("asd;asjd;fj");
  INFO_LOG("asd;asjd;fj");
  DEBUG_LOG("asd;asjd;fj");

  MatrixXd a(3,2);
  a << 12,3,4,5,6,7;
  cout<< endl << a << endl;
  assert(EIGEN3EXT::write("a.b",a));

  Timer timer;
  timer.start();
  MatrixXd b;
  assert(EIGEN3EXT::load("a.b",b));
  cout<< endl << b << endl;
  timer.stop();

  return 0;
}
