#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <PartialConstraints.h>
using namespace Eigen;
using namespace UTILITY;

BOOST_AUTO_TEST_SUITE(PartialConstraintsTest)

BOOST_AUTO_TEST_CASE(testIO){
  
  PartialConstraintsSet parcons;
  TEST_ASSERT(parcons.load(string(TEST_DATA_DIR)+"PartialConTest.txt"));
  ASSERT_EQ( parcons.getPartialConSet().size(),1);
  pPartialConstraints_const p = parcons.getPartialCon(0);
  TEST_ASSERT( p!= NULL);
  if(p){
	ASSERT_EQ( p->numConNodes(),61 );
	ASSERT_EQ( *(p->getConNodesSet()[0].begin()), 170 );
  }
}

BOOST_AUTO_TEST_SUITE_END()
