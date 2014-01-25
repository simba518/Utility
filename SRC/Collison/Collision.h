#ifndef _COLLISION_H_
#define _COLLISION_H_

#include <CollisionPrimaries.h>
#include <assertext.h>

namespace UTILITY{

#define COLLISION_OBJ_NOTFOUND -1
  
  /**
   * @class Collision Interfaces for collision detection and response of
   * deformable objects. The input are nodes and faces of triangle meshes.
   * 
   */
  class Collision{
	
  public:
	virtual int addObject(const void*object_id,const VectorV3&nodes,const VectorV3i&face){

	  assert ( !containObject(object_id) );
	  forces.push_back ( pCollisionForce(new CollisionForce() ) );
	  objects.push_back( object_id );
	  return (int)objects.size();
	}
	template<class VECTOR_D>
	int addObjectV(const void*object_id,const VECTOR_D&nodes,const VectorV3i&face){
	  assert_eq(nodes.size()%3,0);
	  VectorV3 n(nodes.size()/3);
	  for (size_t i = 0; i < n.size(); ++i){
		n[i][0] = nodes[i*3+0];
		n[i][1] = nodes[i*3+1];
		n[i][2] = nodes[i*3+2];
	  }
	  return this->addObject(object_id, n, face);
	}
	template<class VECTOR_D, class VECTOR_I>
	int addObjectV(const void*object_id,const VECTOR_D&nodes,const VECTOR_I&face){
	  assert_eq(face.size()%3,0);
	  VectorV3i f(face.size()/3);
	  for (size_t i = 0; i < f.size(); ++i){
		f[i][0] = face[i*3+0];
		f[i][1] = face[i*3+1];
		f[i][2] = face[i*3+2];
	  }
	  return this->addObjectV(object_id, nodes, f);
	}

	virtual bool prepare() = 0;
	virtual bool update(const void *object_id, const VectorV3& nodes) = 0;
	template<class VECTOR_D>
	bool updateV(const void *object_id, const VECTOR_D& nodes){
	  assert_eq(nodes.size()%3,0);
	  VectorV3 n(nodes.size()/3);
	  for (size_t i = 0; i < n.size(); ++i){
		n[i][0] = nodes[i*3+0];
		n[i][1] = nodes[i*3+1];
		n[i][2] = nodes[i*3+2];
	  }
	  this->update(object_id, n);
	}

	virtual void checkCollision() = 0;
	virtual void computeForce() = 0;
	pCollisionForce_const getForce(const void *object_id)const{
	  const int id = getObjectId(object_id);
	  return getForce(id);
	}
	pCollisionForce_const getForce(const int object_id)const{
	  pCollisionForce_const p_force;
	  if (object_id >= 0&& object_id < (int)forces.size()){
		p_force = forces[object_id];
	  }
	  return p_force;
	}
	bool addCollisionForce(const void *object_id, VectorX &f_ext)const{

	  bool succ = false;
	  pCollisionForce_const fc = getForce(object_id);
	  if (fc){
		succ = true;
		assert_eq (fc->vertex.size(), fc->forces.size());
		for (int i = 0; i < (int)fc->vertex.size(); ++i){
		  const int v = fc->vertex[i];
		  const Vector3d &f = fc->forces[i];
		  assert_in (v*3,0,f_ext.size()-2);
		  f_ext[3*v + 0] += f[0];
		  f_ext[3*v + 1] += f[1];
		  f_ext[3*v + 2] += f[2];
		}
	  }
	  return succ;
	}
	int getObjectId(const void *object_id)const{

	  int id = COLLISION_OBJ_NOTFOUND;
	  for (int i = 0; i < (int)objects.size(); ++i){
		if (objects[i] == object_id){
		  id = i;
		  break;
		}
	  }
	  return id;
	}
	virtual void removeAllCollisionObjects(){
	  forces.clear();
	  objects.clear();
	}
	bool containObject(const void *object_id)const{
	  return containObject( getObjectId(object_id) );
	}
	bool containObject(const int object_id)const{
	  return ( object_id >=0 && object_id < (int)objects.size() );
	}

  protected:
	vector<pCollisionForce> forces;
	vector<const void *> objects;
  };
  
  typedef boost::shared_ptr<Collision> pCollision;  

}//end of namespace

#endif /*_COLLISION_H_*/
