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
#include <MatrixTools.h>
#include <SparseMatrixTools.h>
#include <eigen3/Eigen/Dense>
#include <SparseGenEigenSolver.h>
#include <CASADITools.h>
using namespace UTILITY;
using namespace std;
using namespace Eigen;

int main(int argc, char *argv[]){

  JsonFilePaser file;
  int a;
  std::vector<std::string> av;
  file.read("hh",av);

  return 0;
}
