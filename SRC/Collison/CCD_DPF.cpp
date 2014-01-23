#include <ccdAPI.h>
#include "CCD_DPF.h"
using namespace DEF_COLLISION;

pCCD_DPF CCD_DPF::p_instance;

void EECallback(unsigned int e1_v1, unsigned int e1_v2,
				unsigned int e2_v1, unsigned int e2_v2, float t) {

  // printf("EE result: e1(%d, %d), e2(%d, %d) @ t=%f\n", e1_v1, e1_v2, e2_v1, e2_v2, t);
  EECollisionSet &ees = CCD_DPF::getInstance()->getEECollisionSet();
  ees.insert ( EECollision (e1_v1, e1_v2, e2_v1, e2_v2, t) );
}

void VFCallback(unsigned int vid, unsigned int fid, float t) {

  // printf("VF result: v=%d, f=%d @ t=%f\n", vid, fid, t);
  VFCollisionSet &vfs = CCD_DPF::getInstance()->getVFCollisionSet();
  vfs.insert ( VFCollision (vid,fid,t) );
}

int CCD_DPF::addObject(const void *object_id, const VectorV3& nodes, const VectorV3i &faces){

  Collision::addObject(object_id, nodes, faces);

  int start = -1;
  start = _verts.size();
  _start_index.push_back(start);

  _verts.reserve(_verts.size()+nodes.size());
  _faces.reserve(_faces.size()+faces.size());
  for (size_t i = 0; i < nodes.size(); ++i){
  	_verts.push_back( SELF_CCD::vec3f(nodes[i].x(),nodes[i].y(),nodes[i].z()));
  }
  for (size_t i = 0; i < faces.size(); ++i){
  	_faces.push_back( SELF_CCD::tri3f(start + faces[i].x(),start + faces[i].y(),start + faces[i].z()));
  }

  return start;
}

bool CCD_DPF::prepare(){

  if (_verts.size() > 0){

	ccdInitModel(_verts,_faces);
	ccdSetEECallback(EECallback);
	ccdSetVFCallback(VFCallback);
  }
  return true;
}

bool CCD_DPF::update(const void *object_id, const VectorV3& nodes){

  bool succ = false;
  const int id = Collision::getObjectId(object_id);
  assert (id != COLLISION_OBJ_NOTFOUND);

  if ( id != COLLISION_OBJ_NOTFOUND ){
	
  	const int i0 = _start_index[id];
  	for (size_t i = 0; i < nodes.size(); ++i) {
  	  _verts[i0+i].set_value( nodes[i].x(), nodes[i].y(), nodes[i].z() );
  	}
  	ccdUpdateVtxs(_verts);
  	succ = true;
  }
  return succ;
}

void CCD_DPF::computeForce(){

  EECollisionSet &ees = CCD_DPF::getInstance()->getEECollisionSet();
  VFCollisionSet &vfs = CCD_DPF::getInstance()->getVFCollisionSet();
  cout << "ees: " << ees.size() << endl;
  cout << "vfs: " << vfs.size() << endl;
  
  /// @todo compute collision forces
}

void CCD_DPF::removeAllCollisionObjects(){

  Collision::removeAllCollisionObjects();
  _start_index.clear();
  _verts.clear();
  _faces.clear();
  _EECollisionSet.clear();
  _VFCollisionSet.clear();
  ccdQuitModel();
}

void CCD_DPF::checkCollision(){

  ccdChecking(true);
}
