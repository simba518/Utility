#ifndef _CCD_DPF_H_
#define _CCD_DPF_H_

#include <boost/shared_ptr.hpp>
#include <Collision.h>
#include <ccdAPI.h>

namespace UTILITY{

  class CCD_DPF;
  typedef boost::shared_ptr<CCD_DPF> pCCD_DPF;
  
  /**
   * @class CCD_DPF Continuous Collision Detection with Discrete Penalty Forces.
   * 
   */
  class CCD_DPF: public Collision{
	
  public:
	static pCCD_DPF getInstance() {
	  
	  if(p_instance == NULL){
		p_instance = pCCD_DPF(new CCD_DPF());
	  }
	  return p_instance;
	}
	void setReponseScalor(const double rs){
	  assert_gt(rs,0.0f);
	  response_scalor = rs;
	}
	int addObject(const void *object_id, const VectorV3& nodes, const VectorV3i &faces);
	bool prepare();
	bool update(const void *object_id, const VectorV3& nodes);
	void checkCollision();
	void computeForce();
	void removeAllCollisionObjects();
	~CCD_DPF(){
	  removeAllCollisionObjects();
	}

	const vector<int> &startIndex()const{
	  return _start_index;
	}
	const SELF_CCD::vec3f_list &verts()const{
	  return _verts;
	}
	const SELF_CCD::tri_list &faces()const{
	  return _faces;
	}
	const EECollisionSet &getEECollisionSet()const{
	  return _EECollisionSet;
	}
	const VFCollisionSet &getVFCollisionSet()const{
	  return _VFCollisionSet;
	}
	EECollisionSet &getEECollisionSet(){
	  return _EECollisionSet;
	}
	VFCollisionSet &getVFCollisionSet(){
	  return _VFCollisionSet;
	}

  protected:
	CCD_DPF(){
	  response_scalor = 1.0f;
	}
	void convertToObejctId(const int v, int &obj_v, int &obj_id)const;

  private:
	double response_scalor;
	vector<int> _start_index;
	SELF_CCD::vec3f_list _verts;
	SELF_CCD::tri_list _faces;
	EECollisionSet _EECollisionSet;
	VFCollisionSet _VFCollisionSet;
	static pCCD_DPF p_instance;

  };
    
}//end of namespace

#endif /*_CCD_DPF_H_*/
