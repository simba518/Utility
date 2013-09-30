#ifndef _TETMESHEMBEDING_H_
#define _TETMESHEMBEDING_H_

#include <boost/shared_ptr.hpp>
#include <TetMesh.h>

namespace UTILITY{

  struct VertexMap{

	size_t _tid; // tet id
	Vector4d _bcCoord; // barycenter coordinates.
	// The cloest point on the tet with respect to the corresponding vertex on
	// the mesh. If the vertex is in the tetrahedron, _close is it self,
	// otherwise it is on the surface of the tetrahedron.
	Vector3d _close; 
  };

  // attach each cubes with a set of tets interact with it, use it to speed up
  // the search process.
  typedef boost::unordered_map<HashedId,std::vector<size_t>,HashedIdHash> VertexHash;
  
  /**
   * @class TetMeshEmbeding obj mesh embedded by a tet mesh.
   * 
   */
  class TetMeshEmbeding{

  public:
	TetMeshEmbeding(pTetMesh_const vol_mesh):_vol_mesh(vol_mesh){}
	void embedMesh(const VVec3d &verts);
	void generateNewVertexList(const VectorXd& nodes,VVec3d& verts);
	void interpolateForces(const vector<int>& node_ids,const VVec3d& forces, VectorXd& rlst_forces);
	const std::vector<VertexMap> &mapping()const{
	  return _mapping;
	}

	// io
	bool write(std::ostream& os)const{
	  return false;
	}
	bool read(std::istream& is){
	  return false;
	}

  protected:
	void generateHash(VertexHash& hash,double& w);
	void find(Vector3d& cp,size_t& tid,Vector4d& bary,const size_t& rad,
			  const Vector3d& pos,const VertexHash& hash,const double& w)const;

  private:
	pTetMesh_const _vol_mesh;
	std::vector<VertexMap> _mapping;
  };
  
  typedef boost::shared_ptr<TetMeshEmbeding> pTetMeshEmbeding;

}//end of namespace

#endif /* _TETMESHEMBEDING_H_ */
