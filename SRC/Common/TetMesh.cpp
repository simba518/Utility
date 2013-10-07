#include <float.h>
#include <boost/filesystem/fstream.hpp>
#include <boost/foreach.hpp>
#include <TetMesh.h>
#include <AuxTools.h>
#include <Log.h>
using namespace UTILITY;

/*each face belongs to two elements*/
#define ADD_NEIGH_FACE(a,b,c)											\
  if((iter=_faceId.find(HashedId(tet[a],tet[b],tet[c],0))) == _faceId.end()) \
	_faceId[HashedId(tet[a],tet[b],tet[c],0)]=std::pair<size_t,size_t>(i,-1); \
  else iter->second.second=i;

#define SET_NEIGH_NODE(a,b)						\
  _nodeNeighNode[tet[a]].insert(tet[b]);		\
  _nodeNeighNode[tet[b]].insert(tet[a]);

#define SET_NEIGH(a,b,c,m)												\
  iter=_faceId.find(HashedId(tet[a],tet[b],tet[c],0));					\
  if(iter->second.second == -1)											\
	fn[m]=-1;															\
  else fn[m]=iter->second.first == i ? iter->second.second : iter->second.first;

void TetMesh::reset(const VVec3d& nodes, const VVec4i& tets){

  FaceId::iterator iter;
  _nodes=nodes;
  _tets=tets;
		
  _nodeNeighNode.resize(_nodes.size());
  for(size_t i=0;i<(size_t)_tets.size();i++) {

  	Vector4i& tet=_tets[i];
  	if(tetrahedron(_nodes[tet.x()],_nodes[tet.y()],_nodes[tet.z()],_nodes[tet.w()])._swap)
  	  std::swap(tet.z(),tet.w());

  	ADD_NEIGH_FACE(0,1,2);
  	ADD_NEIGH_FACE(0,2,3);
  	ADD_NEIGH_FACE(0,3,1);
  	ADD_NEIGH_FACE(1,2,3);

  	SET_NEIGH_NODE(0,1);
  	SET_NEIGH_NODE(0,2);
  	SET_NEIGH_NODE(0,3);
  	SET_NEIGH_NODE(1,2);
  	SET_NEIGH_NODE(2,3);
  	SET_NEIGH_NODE(3,1);
  }

  _surface.clear();
  for(iter=_faceId.begin();iter!=_faceId.end();iter++){
    if (iter->second.second < 0){
  	  _surface.push_back(Vector3i( iter->first._id[0], iter->first._id[1], iter->first._id[2] ));
  	}
  }

  _faceNeighTet.resize(_tets.size());
  for(size_t i=0;i<(size_t)_tets.size();i++) {

  	const Vector4i& tet=_tets[i];
  	Vector4i& fn=_faceNeighTet[i];
  	SET_NEIGH(1,2,3,0);
  	SET_NEIGH(0,2,3,1);
  	SET_NEIGH(0,1,3,2);
  	SET_NEIGH(0,1,2,3);
  }
  _mtl.reset(this->tets().size());
}

int TetMesh::getContainingElement(const Vector3d &pos)const{

  // linear scan
  const int numElements = _tets.size();
  for(int element=0; element < numElements; element++){
	const tetrahedron t = getTet(element);
	if (t.isInside(pos)) return element;
  }
  return -1;
}

int TetMesh::getClosestElement(const Vector3d &pos)const{

  const int numElements = _tets.size();
  double closestDist = DBL_MAX;
  int closestElement = 0;
  for(int element=0; element < numElements; element++){
	const tetrahedron t = getTet(element);
	const Vector3d center = t.center();
	const double dist = (pos-center).norm();
	if (dist < closestDist){
	  closestDist = dist;
	  closestElement = element;
	}
  }
  return closestElement;
}

int TetMesh::buildInterpWeights(const VectorXd &vertices,vector<int> &nodes,
								VectorXd &weights,const double zeroThreshold)const{
  TRACE_FUN();

  const int numElementVertices = 4;
  const int numTargetLocations = vertices.size()/3;
  nodes.resize( numElementVertices*numTargetLocations );
  weights.resize( numElementVertices*numTargetLocations );

  Vector4d baryWeights;
  int numExternalVertices = 0;

  INFO_LOG("Wait until reach " << numTargetLocations);
  for (int i=0; i < numTargetLocations; i++){

    if (i%100 == 0) { printf("%d ", i); fflush(NULL);}
    const Vector3d &pos = vertices.segment(i*3,3);
    int element = getContainingElement(pos);
    if (element < 0) {
      element = getClosestElement(pos);
      numExternalVertices++;
    }
	const tetrahedron t = getTet(element);
	baryWeights = t.bary(pos);
    if (zeroThreshold > 0) {
      // check whether vertex is close enough to the mesh
      double minDistance = DBL_MAX;
      int assignedZero = 0;
      for(int ii=0; ii< numElementVertices;ii++) {
		const double dist = (node(element,ii)-pos).norm();
		minDistance = minDistance>dist ? dist:minDistance;
      }
      if (minDistance > zeroThreshold) {
		baryWeights.setZero();
        assignedZero++;
        continue;
      }
    }
    for(int ii=0; ii<numElementVertices; ii++){
      nodes[numElementVertices*i+ii] = _tets[element][ii];
      weights[numElementVertices*i+ii] = baryWeights[ii];
    }
  }
  return numExternalVertices;
}

void TetMesh::interpolate(const vector<int> &tetNodes,const VectorXd &weights,
						  const VectorXd& u,VectorXd& uTarget){

  assert_eq(tetNodes.size(),weights.size());
  assert_eq(tetNodes.size()%4,0);
  const int numTargetLocations = tetNodes.size()/4;
  const int numElementVertices = 4;
  uTarget.resize(numTargetLocations*3);
  Vector3d defo;
  for (int i=0; i < numTargetLocations; i++) {
	defo.setZero();
	for (int j=0; j < numElementVertices; j++) {
	  const int k = tetNodes[numElementVertices*i+j];
	  assert_in(k*3,0,u.size()-2);
	  defo += weights[numElementVertices*i+j]*u.segment(k*3,3);
	}
	uTarget.segment(i*3,3) = defo;
  }
}

bool TetMesh::load(const std::string& filename){

  VVec3d nodes;
  VVec4i tets;
  string line;
  bool succ = false;

  // load nodes
  INFILE(is,filename.c_str());
  while(getline(is,line) && (succ=(line.find("NODE")==string::npos))){};
  Vector3d v;
  char tc;
  while (is.good() && !is.eof()){
	is >> line;
	if (line.find("ELEMENT") != string::npos){
	  succ = true;
	  getline(is,line);
	  break;
	}
	is >> v[0] >> tc >> v[1]>> tc >> v[2];
	nodes.push_back(v);
  }

  // load tets
  Vector4i t;
  const Vector4i t1 = Vector4i::Ones(4);
  succ = false;
  while (is.good() && !is.eof()){
	is >> line;
	if (line.find("ELSET") != string::npos){
	  succ = true;
	  break;
	}
	is >> t[0]>> tc >> t[1]>> tc >> t[2]>> tc >> t[3];
	t -= t1;
	assert_in(t[0],0,nodes.size()-1);
	assert_in(t[1],0,nodes.size()-1);
	assert_in(t[2],0,nodes.size()-1);
	assert_in(t[3],0,nodes.size()-1);
	tets.push_back(t);
  }

  is.close();
  reset(nodes,tets);
  return succ;
}

bool TetMesh::write(const std::string& filename)const{

  return false;
}
