#ifndef _OBJMESH_H_
#define _OBJMESH_H_

#include <map>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <Eigen/Dense>
#include <BBox.h>
#include <assertext.h>
using namespace std;

namespace UTILITY{

  typedef struct{
	
	std::string name;

    float ambient[3];
    float diffuse[3];
    float specular[3];
    float transmittance[3];
    float emission[3];
    float shininess;
    float ior;  // index of refraction

	string ambient_texname;
    string diffuse_texname;
    string specular_texname;
    string normal_texname;
    map<string, string> unknown_parameter;

  }ObjMtl;

  class ObjMesh{
	
  public:
	ObjMesh(){}
	template<class VECTOR_D, class VECTOR_I>
	ObjMesh(const VECTOR_D &vertices,const VECTOR_I &faces){
	  assert_eq(vertices.size()%3,0);
	  setVerts(vertices);
	  setFaces(faces);
	}
	
	// set
	void setMtl(const ObjMtl &mtl){
	  _mtl = mtl;
	}
	void setVerts(const Eigen::VectorXd &v){ _verts = v;}
	template<class VECTOR>
	void setVerts(const VECTOR &v){setVec(v,_verts);}
	template<class VECTOR>
	void setFaces(const VECTOR &f){setVec(f,_faces);}
	template<class VECTOR>
	void setVertNormals(const VECTOR &v){setVec(v,_vertNormal);}

	// get
	const ObjMtl &getMtl()const{
	  return _mtl;
	}
	const Eigen::VectorXd &getVerts()const{
	  return _verts;
	}
	const Eigen::VectorXi &getFaces()const{
	  return _faces;
	}
	const Eigen::VectorXd &getVertNormal()const{
	  return _vertNormal;
	}
	const Eigen::Vector3d getVerts(const int i)const{
	  return getSubV3(_verts,i);
	}
	const Eigen::Vector3i getFaces(const int i)const{
	  return getSubV3(_faces,i);
	}
	const Eigen::Vector3d getVertNormal(const int i)const{
	  return getSubV3(_vertNormal,i);
	}
	int getVertsNum()const{
	  return _verts.size()/3;
	}
	int getFacesNum()const{
	  return _faces.size()/3;
	}

	BBox<double,Eigen::Vector3d> getBBox()const{
	  BBox<double,Eigen::Vector3d> box;
	  if(_verts.size()>0){
		assert_eq(_verts.size()%3,0);
		box.reset(&(_verts[0]),_verts.size()/3);
	  }
	  return box;
	}
	Eigen::Vector3d getCenter()const{
	  Eigen::Vector3d center;
	  center.setZero();
	  for (int i = 0; i < getVertsNum(); ++i)  center += getVerts(i);
	  if(getVertsNum()>0)	center = center/getVertsNum();
	  return center;
	}

  protected:
	template<class T> 
	Eigen::Matrix<T,3,1> getSubV3(const Eigen::Matrix<T,-1,1> &V,int i)const{
	  assert_in(i,0,V.size()/3-1);
	  return V.segment(i*3,3);
	}
	template<class VECTOR, class T> 
	void setVec(const VECTOR &inputV, Eigen::Matrix<T,-1,1> &V)const{
	  V.resize(inputV.size());
	  for (int i = 0; i < V.size(); ++i) V[i] = inputV[i];
	}
	
  private:
	Eigen::VectorXd _verts;
	Eigen::VectorXi _faces;
	Eigen::VectorXd _vertNormal;
	ObjMtl _mtl;
  };
  
  typedef boost::shared_ptr<ObjMesh> pObjMesh;
  typedef boost::shared_ptr<const ObjMesh> pObjMesh_const;

}

#endif /*_OBJMESH_H_*/
