#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <MatrixTools.h>
using namespace Eigen;
using namespace EIGEN3EXT;

BOOST_AUTO_TEST_SUITE(MatrixTools)

BOOST_AUTO_TEST_CASE(createTest){

  const double A[] = {1,2,3,4,5,6};

  MatrixXd m;
  createFromRowMajor(m, &A[0], 2, 3);
  MatrixXd m2(2,3);
  m2 << 1,2,3,4,5,6;
  ASSERT_EQ (m,m2);

  createFromColMajor(m, &A[0], 2, 3);
  m2 << 1,3,5,2,4,6;
  ASSERT_EQ (m,m2);
}

BOOST_AUTO_TEST_CASE(SVDTest){
  
  const Matrix3d M = Matrix3d::Random();
  Matrix3d U, Vt, D;
  ModifiedSVD3x3(M,U,Vt,D);
  ASSERT_EQ_TOL (U.determinant(), 1, 1e-10);
  ASSERT_EQ_TOL (Vt.determinant(), 1, 1e-10);
  ASSERT_EQ_SMALL_MAT_TOL ((U.transpose()*U), Matrix3d::Identity(), 1e-10);
  ASSERT_EQ_SMALL_MAT_TOL ((Vt.transpose()*Vt), Matrix3d::Identity(), 1e-10);
  ASSERT_EQ_SMALL_MAT_TOL ((U*D*Vt), M, 1e-10);
}

BOOST_AUTO_TEST_CASE(PolarDecompositionTest){

  const Matrix3d M = Matrix3d::Random();
  Matrix3d R,S;
  ModifiedPD3x3(M,R,S);
  ASSERT_EQ_TOL (R.determinant(), 1, 1e-10);
  ASSERT_EQ_SMALL_MAT_TOL ((R.transpose()*R), Matrix3d::Identity(), 1e-10);
  ASSERT_EQ_SMALL_MAT_TOL (S.transpose(), S, 1e-10);
  ASSERT_EQ_SMALL_MAT_TOL ((R*S), M, 1e-10);
}

BOOST_AUTO_TEST_CASE(convertTest){
  
  MatrixXd U(2,3);
  std::vector<VectorXd> u;
  convert(U,u);
  ASSERT_EQ(u.size(),U.cols());
  ASSERT_EQ(u[0].size(),U.rows());
  ASSERT_EQ(u[1].size(),U.rows());
  ASSERT_EQ(u[2].size(),U.rows());
  
  MatrixXd M;
  convert(u,M);
  ASSERT_EQ(M.rows(),U.rows());
  ASSERT_EQ(M.cols(),U.cols());

  ASSERT_EQ(M,U);
}

BOOST_AUTO_TEST_SUITE_END()
