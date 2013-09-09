#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>

#include <eigen3/Eigen/Dense>
using namespace Eigen;

BOOST_AUTO_TEST_SUITE(FancyShmancyLogic2)

BOOST_AUTO_TEST_CASE(TestingIf2x3equals600)
{
  MatrixXd M(1,2), A(1,2);
  M << 1,2;
  A << 3,4;
  ASSERT_EQ_SMALL_MAT(M,A);
}

BOOST_AUTO_TEST_CASE(TestingIf2x2equals4)
{
}

BOOST_AUTO_TEST_SUITE_END()
