#include <iostream>
#include <MatrixTools.h>
#include <eigen3/Eigen/Dense>
using namespace std;
using namespace Eigen;

int main(int argc, char *argv[]){
  
  MatrixXd M1 = MatrixXd::Random(4,4)+MatrixXd::Identity(4,4)*10.0f;
  const MatrixXd M = M1.transpose() + M1;
  MatrixXd U = MatrixXd::Random(4,2);

  cout << U.transpose()*M*U << endl;
  EIGEN3EXT::MGramSchmidt(M,U);
  cout << U.transpose()*M*U << endl;
  return 0;
}
