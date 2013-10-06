#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <TetMesh.h>
#include <TetMeshEmbeding.h>
using namespace Eigen;
using namespace UTILITY;

/**
 * @class TetMeshFactoryForTest generate tet meshes for test.
 * 
 */
class TetMeshFactoryForTest{
	
public:
  static pTetMesh tet2(){ 

	// create a simple mesh with two tets with five nodes
	VVec3d nodes;
	VVec4i tets;
	nodes.push_back(Vector3d(0,0,1));
	nodes.push_back(Vector3d(1,0,0));
	nodes.push_back(Vector3d(0,0,0));
	nodes.push_back(Vector3d(0,1,0));
	nodes.push_back(Vector3d(0,-5,0));

	tets.push_back(Vector4i(0,1,2,3));
	tets.push_back(Vector4i(0,2,1,4));

	pTetMesh tet_mesh= pTetMesh(new TetMesh(nodes, tets));
	tet_mesh->setSingleMaterial(50.0f,0.1f,0.45f);
	return tet_mesh;
  }

  static pTetMesh tet1(){

	// create a simple mesh with two tets with five nodes
	VVec3d nodes;
	VVec4i tets;
	nodes.push_back(Vector3d(0,0,1));
	nodes.push_back(Vector3d(1,0,0));
	nodes.push_back(Vector3d(0,0,0));
	nodes.push_back(Vector3d(0,1,0));

	tets.push_back(Vector4i(0,1,2,3));

	pTetMesh tet_mesh= pTetMesh(new TetMesh(nodes, tets));
	tet_mesh->setSingleMaterial(50.0f,0.5,0.45);
	return tet_mesh;
  }
};

BOOST_AUTO_TEST_SUITE(TetMeshTest)

BOOST_AUTO_TEST_CASE(initializeTest){

  // create a simple mesh with two tets with five nodes
  pTetMesh tet_mesh = TetMeshFactoryForTest::tet2();
  TEST_ASSERT (tet_mesh != NULL);

  // check nodes, tets
  ASSERT_EQ (tet_mesh->nodes().size(), 5);
  ASSERT_EQ (tet_mesh->tets().size(), 2);
  { // check face ids
  	const FaceId &faceId = tet_mesh->faceId();
  	ASSERT_EQ (faceId.size(), 7);
  	FaceId::const_iterator iter;
  	iter = faceId.find(HashedId(1,2,3,0));
  	TEST_ASSERT ( iter != faceId.end() );
  	ASSERT_EQ (iter->second.first,0);
  	ASSERT_EQ (iter->second.second,-1);

  	iter=faceId.find(HashedId(0,2,3,0));
  	TEST_ASSERT ( iter != faceId.end() );
  	ASSERT_EQ (iter->second.first,0);
  	ASSERT_EQ (iter->second.second,-1);

  	iter=faceId.find(HashedId(0,1,3,0));
  	TEST_ASSERT ( iter != faceId.end() );
  	ASSERT_EQ (iter->second.first,0);
  	ASSERT_EQ (iter->second.second,-1);

  	iter=faceId.find(HashedId(0,1,2,0));
  	TEST_ASSERT ( iter != faceId.end() );
  	ASSERT_EQ (iter->second.first,0);
  	ASSERT_EQ (iter->second.second,1);

  	iter=faceId.find(HashedId(1,2,4,0));
  	TEST_ASSERT ( iter != faceId.end() );
  	ASSERT_EQ (iter->second.first,1);
  	ASSERT_EQ (iter->second.second,-1);

  	iter=faceId.find(HashedId(0,2,4,0));
  	TEST_ASSERT ( iter != faceId.end() );
  	ASSERT_EQ (iter->second.first,1);
  	ASSERT_EQ (iter->second.second,-1);

  	iter=faceId.find(HashedId(0,1,4,0));
  	TEST_ASSERT ( iter != faceId.end() );
  	ASSERT_EQ (iter->second.first,1);
  	ASSERT_EQ (iter->second.second,-1);

  	iter=faceId.find(HashedId(0,1,8,0));
  	TEST_ASSERT ( iter == faceId.end() );
  }
  
  { // check neighors of tets sharing the same faces
  	const VVec4i &faceNeighTet = tet_mesh->faceNeighTet();
  	ASSERT_EQ (faceNeighTet.size(),2);
  	ASSERT_EQ ( faceNeighTet[0], Vector4i(-1,-1,-1,1));
  	ASSERT_EQ ( faceNeighTet[1], Vector4i(-1,-1,-1,0));
  }
  
  { // check neighors of nodes share the same edges.
  	const VectorUseti &nodeNeighNode = tet_mesh->nodeNeighNode();

  	ASSERT_EQ( nodeNeighNode.size(), 5 );
  	ASSERT_EQ( nodeNeighNode[0].size(), 4 );
  	ASSERT_EQ( nodeNeighNode[1].size(), 4 );
  	ASSERT_EQ( nodeNeighNode[2].size(), 4 );
  	ASSERT_EQ( nodeNeighNode[3].size(), 3 );
  	ASSERT_EQ( nodeNeighNode[4].size(), 3 );

  	TEST_ASSERT( nodeNeighNode[0].find(1) != nodeNeighNode[0].end() );
  	TEST_ASSERT( nodeNeighNode[0].find(2) != nodeNeighNode[0].end() );
  	TEST_ASSERT( nodeNeighNode[0].find(3) != nodeNeighNode[0].end() );
  	TEST_ASSERT( nodeNeighNode[0].find(4) != nodeNeighNode[0].end() );
  }

  // test empty mesh -----------------------------------------------------------
  VVec3d nodes;
  VVec4i tets;
  tet_mesh->reset(nodes, tets);
  ASSERT_EQ (tet_mesh->nodes().size(), 0);
}

BOOST_AUTO_TEST_CASE(IOTest){
  
}

BOOST_AUTO_TEST_CASE(testContainingEle){
 
  pTetMesh tet_mesh = TetMeshFactoryForTest::tet2();
  Vector3d pos;

  pos << -1,0,0;
  ASSERT_EQ (tet_mesh->getContainingElement(pos),-1);

  pos << 0,0,0;
  ASSERT_EQ (tet_mesh->getContainingElement(pos),0);

  pos << 0,-1,0;
  ASSERT_EQ (tet_mesh->getContainingElement(pos),1);

  pos << 0.1,-0.1,0.1;
  ASSERT_EQ (tet_mesh->getContainingElement(pos),1);

  pos << 0,-6,0;
  ASSERT_EQ (tet_mesh->getContainingElement(pos),-1);
}

BOOST_AUTO_TEST_CASE(testClosestEle){

  pTetMesh tet_mesh = TetMeshFactoryForTest::tet2();
  Vector3d pos;

  pos << 10,0,0;
  ASSERT_EQ (tet_mesh->getClosestElement(pos),0);

  pos << -1,0,0;
  ASSERT_EQ (tet_mesh->getClosestElement(pos),0);

  pos << 0,0,0;
  ASSERT_EQ (tet_mesh->getClosestElement(pos),0);

  pos << 0,-1,0;
  ASSERT_EQ (tet_mesh->getClosestElement(pos),1);

  pos << 0.1,-0.1,0.1;
  ASSERT_EQ (tet_mesh->getClosestElement(pos),0);

  pos << 0,-6,0;
  ASSERT_EQ (tet_mesh->getClosestElement(pos),1);
}

BOOST_AUTO_TEST_CASE(testInterp){
 
  pTetMesh tet_mesh = TetMeshFactoryForTest::tet2();

  VectorXd vertices;
  vector<int> nodes;
  VectorXd weights; 
  tet_mesh->buildInterpWeights(vertices,nodes,weights);

  vertices.resize(9);
  vertices << 1,2,3,4,5,6,7,8,9;
  tet_mesh->buildInterpWeights(vertices,nodes,weights);

  const VectorXd u = VectorXd::Random(nodes.size()*3);
  VectorXd uTarget = vertices;

  cout << uTarget.transpose() << endl << endl;
  tet_mesh->interpolate(nodes,weights,u,uTarget);
  cout << uTarget.transpose() << endl;

}

BOOST_AUTO_TEST_SUITE_END()
