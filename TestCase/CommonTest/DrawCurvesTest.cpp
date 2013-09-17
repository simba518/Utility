#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <DrawCurves.h>
using namespace Eigen;
using namespace UTILITY;

BOOST_AUTO_TEST_SUITE(DrawCurvesTest)

BOOST_AUTO_TEST_CASE(testfun){

  const VectorXd v = VectorXd::Random(10);
  TEST_ASSERT( PythonScriptDraw2DCurves::write(".tempt.py",v) );
}

BOOST_AUTO_TEST_SUITE_END()
