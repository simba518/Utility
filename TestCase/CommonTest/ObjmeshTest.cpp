#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <Objmesh.h>
#include <ObjFileIO.h>
using namespace Eigen;
using namespace UTILITY;

BOOST_AUTO_TEST_SUITE(ObjmeshTest)

BOOST_AUTO_TEST_CASE(testLoadObjFile){

  const string fname = "./TestCase/TestData/dino.obj";
  Objmesh mesh;
  TEST_ASSERT( load(fname,mesh) );
  ASSERT_EQ (mesh.getVertsNum(),28098);
  ASSERT_EQ (mesh.getFacesNum(),56192);

  Vector3d v0,vT;
  v0 << 1.79623, -4.35718, 0.464406;
  vT << 3.51123,4.55572, -0.326623;
  ASSERT_EQ_SMALL_VEC_TOL (mesh.getVerts(0),v0,3,1e-6);
  ASSERT_EQ_SMALL_VEC_TOL (mesh.getVerts(mesh.getVertsNum()-1),vT,3,1e-6);

  Vector3d vn0,vnT;
  vn0 << 0.0862053, -0.99618,0.0139087;
  vnT << -0.0341522, 0.996979, 0.0697619;
  ASSERT_EQ_SMALL_VEC_TOL (mesh.getVertNormal(0),vn0,3,1e-6);
  ASSERT_EQ_SMALL_VEC_TOL (mesh.getVertNormal(mesh.getVertsNum()-1),vnT,3,1e-6); 

  Vector3i f0,fT;
  f0 << 0,1,2;
  fT << 28094,28090,28091;
  ASSERT_EQ_SMALL_VEC_TOL (mesh.getFaces(0),f0,3,1e-6);
  ASSERT_EQ_SMALL_VEC_TOL (mesh.getFaces(mesh.getFacesNum()-1),fT,3,1e-6);

  Vector3d Kd,Ka,Tf,Ks;
  Kd << 0.00,0.60,0.00;
  Ka << 0.00,0.10,0.00;
  Ks << 0.35,0.35,0.35;
  const double Ni = 1.00;
  const double Ns = 200;
  
  ASSERT_EQ_SMALL_VEC_TOL(mesh.getMtl().diffuse,Kd,3,1e-6);
  ASSERT_EQ_SMALL_VEC_TOL(mesh.getMtl().ambient,Ka,3,1e-6);
  ASSERT_EQ_SMALL_VEC_TOL(mesh.getMtl().specular,Ks,3,1e-6);
  ASSERT_EQ_TOL(mesh.getMtl().shininess,Ns,1e-6);
  ASSERT_EQ_TOL(mesh.getMtl().ior,Ni,1e-6);
}

BOOST_AUTO_TEST_SUITE_END()
