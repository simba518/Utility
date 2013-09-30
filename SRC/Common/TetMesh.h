#ifndef _TETMESH_H_
#define _TETMESH_H_

#include <boost/shared_ptr.hpp>
#include <eigen3/Eigen/Dense>
#include <algorithm>
#include <HashedId.h>
using namespace std;
using namespace Eigen;

namespace UTILITY{
 
  template <typename T>
  struct ElasticMaterial{

	void reset(const int elements_num){
	  _rho.assign(elements_num,0.1f);
	  _E.assign(elements_num,0.2f);
	  _v.assign(elements_num,0.3f);
	}
  
	std::vector<T> _rho; // Density
	std::vector<T> _E; // Young's modulus
	std::vector<T> _v; // Poisson's ratio
  
  };

  template <typename T>
  class tetrahedronTpl{

  public:
	typedef typename Eigen::Matrix<T,2,1> PT2;
	typedef typename Eigen::Matrix<T,3,1> PT;
	typedef typename Eigen::Matrix<T,4,1> PT4;
	typedef typename Eigen::Matrix<T,6,1> PT6;
	typedef typename Eigen::Matrix<T,3,3> MAT3;
  public:
	tetrahedronTpl(){}
	tetrahedronTpl(const PT& a,const PT& b,const PT& c,const PT& d)
	  :_a(a),_b(b),_c(c),_d(d){
	  _swap=false;
	  if(volume() < 0.0f){
		  swap(_c,_d);
		  _swap=true;
		}
	}
	PT4 bary(const PT& pt) const{
	  MAT3 A;
	  A.col(0)=_a-_d;
	  A.col(1)=_b-_d;
	  A.col(2)=_c-_d;
	  PT abc=A.inverse()*(pt-_d);
	  return PT4(abc.x(),abc.y(),abc.z(),1.0f-abc.sum());
	}
	bool isInside(const PT& pt) const{
	  PT4 bc=bary(pt);
	  return bc.x() >= 0 && bc.y() >= 0 && bc.z() >= 0 && bc.w() >= 0;
	}
	T volume() const{
	  return (_b-_a).cross(_c-_a).dot(_d-_a)/6.0f;
	}
	const PT& getNode(const int& i) const{return (&_a)[i];}
  public:
	//data
	PT _a;
	PT _b;
	PT _c;
	PT _d;
	bool _swap;
  };
  typedef tetrahedronTpl<double> tetrahedron;

  typedef std::vector<Vector4i,Eigen::aligned_allocator<Vector4i> > VVec4i;
  typedef std::vector<Vector3i,Eigen::aligned_allocator<Vector3i> > VVec3i;
  typedef std::vector<Vector3d,Eigen::aligned_allocator<Vector3d> > VVec3d;
  typedef boost::unordered_map<HashedId,std::pair<int,int>,HashedIdHash> FaceId;
  typedef std::vector<boost::unordered_set<int> > VectorUseti;
 
  class TetMesh{

  public:
	// init
	TetMesh(){}
	TetMesh(const VVec3d& nodes, const VVec4i& tets){
	  reset(nodes,tets);
	}
	void reset(const VVec3d& nodes, const VVec4i& tets);

	// set 
	void setSingleMaterial(const double&dens,const double&E,const double&v){
	  _mtl._rho.assign(_tets.size(),dens);
	  _mtl._E.assign(_tets.size(),E);
	  _mtl._v.assign(_tets.size(),v);
	}
	void setRestPos(const VVec3d& pos){
	  _nodes = pos;
	}

	// get 
	const VVec3d& nodes() const{return _nodes;}
	const VVec4i& tets() const{return _tets;}
	const ElasticMaterial<double>& material() const{return _mtl;}
	double volume(const int i) const{
	  return tetrahedron(_nodes[_tets[i][0]],
	  					 _nodes[_tets[i][1]],
	  					 _nodes[_tets[i][2]],
	  					 _nodes[_tets[i][3]]).volume();
	}
	const VVec4i &faceNeighTet()const{
	  return _faceNeighTet;
	}
	const VectorUseti &nodeNeighNode()const{
	  return _nodeNeighNode;
	}

	// io
	bool read(const std::string& filename);
	bool write(const std::string& filename)const;

  private:
	VVec3d _nodes;
	VVec4i _tets;
	ElasticMaterial<double> _mtl;

	FaceId _faceId; // 3 vertices of each face, and 2 tets it belongs to.
	VVec4i _faceNeighTet; // 4 neighbor tets(same face) of one tet.
	VectorUseti _nodeNeighNode; // the neighbor nodes of each node.
	VVec3i _surface; // surface's faces
  };

  typedef boost::shared_ptr<TetMesh> pTetMesh;
  typedef boost::shared_ptr<const TetMesh> pTetMesh_const;
}

#endif /* _TETMESH_H_ */
