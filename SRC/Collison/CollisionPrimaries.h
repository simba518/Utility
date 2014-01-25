#ifndef _COLLISIONPRIMARIES_H_
#define _COLLISIONPRIMARIES_H_

#include <vector>
#include <algorithm>
#include <eigen3/Eigen/Dense>
#include <boost/unordered_set.hpp>
#include <boost/shared_ptr.hpp>
using namespace Eigen;
using namespace std;

namespace UTILITY{

  typedef std::vector<Vector3d,Eigen::aligned_allocator<Vector3d> > VectorV3;
  typedef Eigen::Matrix<double,-1, 1> VectorX;
  typedef std::vector<Vector3i,Eigen::aligned_allocator<Vector3i> > VectorV3i;

  // edge-edge collision
  struct EECollision{

	EECollision (){
	  _e1v1 = -1;
	  _e1v2 = -1;
	  _e2v1 = -1;
	  _e2v2 = -1;
	  _t = -1;
	}

	EECollision (unsigned int e1v1, unsigned int e1v2,
				 unsigned int e2v1, unsigned int e2v2){

	  _e1v1 = e1v1;
	  _e1v2 = e1v2;
	  _e2v1 = e2v1;
	  _e2v2 = e2v2;
	  _t = -1.0f;
	  this->sort (_e1v1, _e1v2, _e2v1, _e2v2);
	}

	EECollision (unsigned int e1v1, unsigned int e1v2,
				 unsigned int e2v1, unsigned int e2v2, float t){

	  _e1v1 = e1v1;
	  _e1v2 = e1v2;
	  _e2v1 = e2v1;
	  _e2v2 = e2v2;
	  _t = t;
	  this->sort (_e1v1, _e1v2, _e2v1, _e2v2);
	}

	friend bool operator==(const EECollision& left, const EECollision& right){

	  const bool e1 = (left.e1v1() == right.e1v1() && left.e1v2() == right.e1v2());
	  const bool e2 = (left.e2v1() == right.e2v1() && left.e2v2() == right.e2v2());
	  return e1 && e2;
	}

	unsigned int e1v1()const{
	  return _e1v1;
	}

	unsigned int e1v2()const{
	  return _e1v2;
	}

	unsigned int e2v1()const{
	  return _e2v1;
	}

	unsigned int e2v2()const{
	  return _e2v2;
	}

	float t()const{
	  return _t;
	}

	static void sort(unsigned int &e1v1, unsigned int &e1v2,
					 unsigned int &e2v1, unsigned int &e2v2){
	  if (e1v1 > e1v2)
		std::swap(e1v1, e1v2);
	  if (e2v1 > e2v2)
		std::swap(e2v1, e2v2);
	  if (e1v2 > e2v2){
		std::swap(e1v1, e2v1);
		std::swap(e1v2, e2v2);
	  }
	}

  private:
	unsigned int _e1v1;
	unsigned int _e1v2;
	unsigned int _e2v1;
	unsigned int _e2v2;
	float _t;
  };
  inline size_t hash_value(const EECollision ee){

	std::size_t seed = 0;
	boost::hash_combine(seed, ee.e1v1());
	boost::hash_combine(seed, ee.e1v2());
	boost::hash_combine(seed, ee.e2v1());
	boost::hash_combine(seed, ee.e2v2());
	return seed;
  }
  typedef boost::unordered_set<EECollision> EECollisionSet;

  // vertex-face collision
  struct VFCollision{

	VFCollision (){
	  _vid = -1;
	  _fid = -1;
	  _t = -1;
	}

	VFCollision (unsigned int vid, unsigned int fid){
	  _vid = vid;
	  _fid = fid;
	  _t = -1;
	}

	VFCollision (unsigned int vid, unsigned int fid, float t){
	  _vid = vid;
	  _fid = fid;
	  _t = t;
	}

	friend bool operator==(const VFCollision& left, const VFCollision& right){
	  return left.vid() == right.vid() && left.fid() == right.fid();
	}

	unsigned int vid()const{
	  return _vid;
	}

	unsigned int fid()const{
	  return _fid;
	}

	float t()const{
	  return _t;
	}

  private:
	unsigned int _vid;
	unsigned int _fid;
	float _t;
  };
  inline size_t hash_value(const VFCollision vf){

	std::size_t seed = 0;
	boost::hash_combine(seed, vf.fid());
	boost::hash_combine(seed, vf.vid());
	return seed;
  }
  typedef boost::unordered_set<VFCollision> VFCollisionSet;

  struct CollisionForce{

	vector<int> vertex;
	VectorV3 forces;
  };
  typedef boost::shared_ptr<CollisionForce>  pCollisionForce;
  typedef boost::shared_ptr<const CollisionForce>  pCollisionForce_const;

}//end of namespace

#endif /* _COLLISIONPRIMARIES_H_ */
