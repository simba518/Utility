#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <SparseGenEigenSolver.h>
#include <iostream>
using namespace std;
using namespace Eigen;
using namespace EIGEN3EXT;

BOOST_AUTO_TEST_SUITE(SparseGenEigenSolverTest)

/**
 * K:
 * 1.0  0.1  0.0 
 * 0.1  2.0  0.0
 * 0.0  0.0  3.0
 *
 * M:
 * 2.3  0.0  0.0 
 * 0.0  1.2  0.0
 * 0.0  0.0  3.9
 */
BOOST_AUTO_TEST_CASE(TestEigenSolver){
  
  const int n = 3;
  Eigen::SparseMatrix<double> K(n,n),M(n,n);
  K.insert (0,0) = 1.0f;
  K.insert (1,0) = 0.10f;
  K.insert (1,1) = 2.0f;
  K.insert (2,2) = 3.0f;

  M.insert (0,0) = 2.3f;
  M.insert (1,1) = 1.2f;
  M.insert (2,2) = 3.9f;

  MatrixXd eig_vec;
  VectorXd eig_val;
  const int max_eig_num = 2;
  K.makeCompressed();
  M.makeCompressed();

  TEST_ASSERT( EigenSparseGenEigenSolver<double>::solve(K,M,eig_vec,eig_val,max_eig_num) );

  const MatrixXd Km = K;
  const MatrixXd Mm = M;
  Eigen::GeneralizedSelfAdjointEigenSolver<MatrixXd > es(Km,Mm);
  const MatrixXd eigen_m = es.eigenvectors().topLeftCorner(n,max_eig_num);
  
  for (int i = 0; i < max_eig_num; ++i){
	
  	TEST_ASSERT( abs(es.eigenvalues()[i] - eig_val[i]) < 1e-8 );
	TEST_ASSERT( ((eigen_m.col(i)-eig_vec.col(i)).norm() < 1e-8) ||  
					((eigen_m.col(i)+eig_vec.col(i)).norm() < 1e-8) );
  }

}

BOOST_AUTO_TEST_SUITE_END()
