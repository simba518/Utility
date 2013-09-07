#define LOG_WARN
#define LOG_ERROR
#define LOG_DEBUG
#define LOG_INFO
#define LOG_TRACE
// #include <Timer.h>
#include <IO.h>
#include <assertext.h>
#include <iostream>
#include <Log.h>
using namespace UTILITY;

int main(int argc, char *argv[]){

  TRACE_FUN();
  ERROR_LOG("asd;asjd;fj");
  ERROR_LOG_COND("asd;asjd;fj",false);
  WARN_LOG("asd;asjd;fj",true);
  DEBUG_LOG("asd;asjd;fj");
  INFO_LOG("asd;asjd;fj");
  DEBUG_LOG("asd;asjd;fj");

  // Timer timer;
  // timer.start();
  // timer.stop();

  // std::cout<<std::setprecision(12) << "j;askdj" << std::endl;
  
  return 0;
}
