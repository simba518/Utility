#include <boost/filesystem/fstream.hpp>
#include <boost/foreach.hpp>
#include <TetMesh.h>
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

bool TetMesh::read(const std::string& filename){

  return false;
}

bool TetMesh::write(const std::string& filename)const{

  return false;
}
