#ifndef _HASHEDID_H_
#define _HASHEDID_H_

#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <eigen3/Eigen/Dense>

namespace UTILITY{
  
  struct HashedId{
  public:
	HashedId(){
	  _id[0]=_id[1]=_id[2]=_id[3]=-1;
	  _no=-1;
	}
	HashedId(int a,int b,const int& no){
	  if(a > b)
		std::swap(a,b);
	  _id[0]=-1;
	  _id[1]=-1;
	  _id[2]=a;
	  _id[3]=b;
	  _no=no;
	}
	HashedId(int a,int b,int c,const int& no){

	  if(a > b)std::swap(a,b);
	  if(b > c)std::swap(b,c);
	  if(a > b)std::swap(a,b);

	  _id[0]=a;
	  _id[1]=b;
	  _id[2]=c;
	  _id[3]=-1;

	  _no=no;
	}
	HashedId(const Eigen::Vector3i& pos,const int& edgeLen,const int& no){

	  _id[0]=edgeLen;
	  _id[1]=pos.x();
	  _id[2]=pos.y();
	  _id[3]=pos.z();
	  _no=no;
	}
	bool operator==(const HashedId& other) const{

	  return _id[0] == other._id[0] &&
		_id[1] == other._id[1] &&
		_id[2] == other._id[2] &&
		_id[3] == other._id[3];
	}
	int _id[4];
	int _no;
  };

  struct HashedIdHash: public boost::hash<HashedId>{

	std::size_t operator()(const HashedId& i) const{
	  const boost::hash<int> h;
	  return h(i._id[0])+h(i._id[1])+h(i._id[2])+h(i._id[3]);
	}
  };
}

#endif /* _HASHEDID_H_ */

