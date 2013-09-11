#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <CASADITools.h>
using namespace Eigen;
using namespace CASADI;

BOOST_AUTO_TEST_SUITE(CASADITools)

BOOST_AUTO_TEST_CASE(convert2VecTest){

  VectorXd v(3);
  v << 1,2,5;
  const CasADi::SXMatrix SV = convert(v);
  ASSERT_EQ(SV.size1(),v.size());
  ASSERT_EQ(SV.size2(),1);
  const VectorXd a = convert2Vec<double>(SV);
  ASSERT_EQ(a,v);
}

BOOST_AUTO_TEST_SUITE_END()
