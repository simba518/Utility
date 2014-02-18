#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <ReducedSimCubature.h>
using namespace Eigen;
using namespace SIMULATOR;

BOOST_AUTO_TEST_SUITE(ReducedSimCubatureTest)

BOOST_AUTO_TEST_CASE(tetRun){

  MatrixXd B;
  pTetMesh_const tet_mesh;
  ReducedSimCubature cuba(B, tet_mesh);

  MatrixXd trainingReduedDisp, trainingFullForces;
  cuba.run(trainingFullForces, trainingFullForces, 0.01f, 100, 10, 10, 10);
}

BOOST_AUTO_TEST_SUITE_END()
