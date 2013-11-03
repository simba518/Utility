#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <MatrixIO.h>
using namespace EIGEN3EXT;

BOOST_AUTO_TEST_SUITE(MatrixIO)

BOOST_AUTO_TEST_CASE(vectorIO){
  Eigen::VectorXd v1(3),v2;
  v1 << 1,2,3;
  const std::string fname = std::string(TEST_DATA_DIR)+"temptV.b";
  TEST_ASSERT(write(fname,v1));
  TEST_ASSERT(load(fname,v2));
  ASSERT_EQ_SMALL_VEC(v1,v2,v1.size());
}

BOOST_AUTO_TEST_CASE(matrixIO){
  Eigen::MatrixXd m1(3,2);
  Eigen::MatrixXd m2 = Eigen::MatrixXd::Random(3,2);
  ASSERT_NE(m1,m2);
  const std::string fname = std::string(TEST_DATA_DIR)+"temptM.b";
  TEST_ASSERT(write(fname,m1,UTILITY::BINARY));
  TEST_ASSERT(load(fname,m2,UTILITY::BINARY));
  ASSERT_EQ_SMALL_MAT(m1,m2);
}

BOOST_AUTO_TEST_SUITE_END()
