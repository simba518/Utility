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
  else fn[m]=iter->second.first == (int)i ? iter->second.second : iter->second.first;

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
	  Vector3i face(iter->first._id[0],iter->first._id[1],iter->first._id[2] );
	  const Vector3d v1 = _nodes[face[1]]-_nodes[face[0]];
	  const Vector3d v2 = _nodes[face[2]]-_nodes[face[0]];
	  const Vector3d n = v1.cross(v2);
	  const Vector4i &t = _tets[iter->second.first];
	  tetrahedronTpl<double> oneTet(_nodes[t[0]],_nodes[t[1]],_nodes[t[2]],_nodes[t[3]]);
	  const Vector3d vc = oneTet.center()-_nodes[face[0]];
	  if(n.dot(vc) > 0){
		const int f = face[0];
		face[0] = face[1];
		face[1] = f;
	  }
  	  _surface.push_back(face);
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

  computeSurfaceNormal();
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

  assert_eq((int)tetNodes.size(),weights.size());
  assert_eq(tetNodes.size()%4,0);
  const int numElementVertices = 4;
  const int numTargetLocations = tetNodes.size()/numElementVertices;
  uTarget.resize(numTargetLocations*3);
  Vector3d defo;
  for (int i=0; i < numTargetLocations; i++) {
	defo.setZero();
	for (int j=0; j < numElementVertices; j++) {
	  const int k = tetNodes[numElementVertices*i+j];
	  assert_in(k*3,0,u.size()-2);
	  defo += weights[numElementVertices*i+j]*u.segment<3>(k*3);
	}
	uTarget.segment<3>(i*3) = defo;
  }
}

void TetMesh::buildInterpMatrix(const vector<int> &tetNodes,const VectorXd &weights,
								const int tetNodesNum, SparseMatrix<double>& A){
  
  const int n = weights.size()/4;
  vector< Triplet<double> > tri;
  tri.reserve(3*tetNodes.size());

  for (int i = 0; i < n; ++i){
	for (int j = i*4; j < i*4+4; ++j){
	  tri.push_back(Triplet<double>(i*3+0, tetNodes[j]*3+0,weights[j] ));
	  tri.push_back(Triplet<double>(i*3+1, tetNodes[j]*3+1,weights[j] ));
	  tri.push_back(Triplet<double>(i*3+2, tetNodes[j]*3+2,weights[j] ));
	}
  }

  A.resize(n*3,tetNodesNum*3);
  A.setFromTriplets(tri.begin(),tri.end());
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
	assert_in(t[0],0,(int)nodes.size()-1);
	assert_in(t[1],0,(int)nodes.size()-1);
	assert_in(t[2],0,(int)nodes.size()-1);
	assert_in(t[3],0,(int)nodes.size()-1);
	tets.push_back(t);
  }

  is.close();
  reset(nodes,tets);
  return succ;
}

bool TetMesh::write(const std::string& filename)const{

  // load nodes
  OUTFILE(outf,filename.c_str());
  if(!outf.is_open()){
	ERROR_LOG("failed to open file " << filename << " for write!");
	return false;
  }
  outf << "*NODE\n";
  for (size_t i = 0; i < _nodes.size(); ++i)
    outf<<i+1<<", "<<_nodes[i][0]<<", "<<_nodes[i][1]<<", "<< _nodes[i][2] << endl;
  outf << "*ELEMENT, TYPE=C3D4\n";
  for (size_t i = 0; i < _tets.size(); ++i){
    outf<<i+1<<", "<<_tets[i][0]<<", "<<_tets[i][1]<<", "<<_tets[i][2]<<", "<<_tets[i][3]<<"\n";
  }
  outf << "*ELSET,ELSET=EALL,GENERATE\n";
  outf <<"1," << _tets.size() << endl;
  return outf.good();
}

BBoxD TetMesh::getBBox()const{

  BBoxD box;
  VectorXd x;
  nodes(x);
  if(x.size()>0){
	assert_eq(x.size()%3,0);
	box.reset(&(x[0]),x.size()/3);
  }
  return box;
}

bool TetMesh::writeVTK(const std::string& filename,const bool binary)const{

  VTKWriter writer(binary);
  writer.addPoints(_nodes);
  writer.addTets(_tets);
  writer.addTriangles(_surface);
  return writer.write(filename);
}

bool TetMesh::writeVTK(const std::string&filename,const MatrixXd &U,const bool binary)const{

  bool succ = true;
  for (int i = 0; i < U.cols() && succ; ++i){
    const VectorXd &u = U.col(i);
	succ = writeVTK(filename+TOSTR(i)+".vtk",u,binary);
  }
  return succ;
}

void TetMesh::computeSurfaceNormal(){
  
  const VVec3d &p = nodes();
  const VVec3i &f = surface();
  _normal.resize(f.size());
  for (int i = 0; i < _normal.size(); ++i){
	assert_in(f[i][0],0,p.size()-1);
	assert_in(f[i][1],0,p.size()-1);
	assert_in(f[i][2],0,p.size()-1);
	const Vector3d &v0 = p[f[i][0]];
	const Vector3d &v1 = p[f[i][1]];
	const Vector3d &v2 = p[f[i][2]];
    _normal[i] = (v1-v0).cross(v2-v0).normalized();
  }
}
