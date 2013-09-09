#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>

BOOST_AUTO_TEST_SUITE(FancyShmancyLogic)

BOOST_AUTO_TEST_CASE(TestingIf2x3equals600)
{
  double a = 1,b = 1;
  ASSERT_NE(a,b);
}

BOOST_AUTO_TEST_CASE(TestingIf2x2equals433)
{
}

BOOST_AUTO_TEST_SUITE_END()
