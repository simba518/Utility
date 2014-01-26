#include <ccdAPI.h>
#include "CCD_DPF.h"
#include <boost/foreach.hpp>
using namespace UTILITY;

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

int CCD_DPF::addObject(const void *object_id,const VectorV3&nodes,const VectorV3i&faces){

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

  const EECollisionSet &ees = CCD_DPF::getInstance()->getEECollisionSet();
  const VFCollisionSet &vfs = CCD_DPF::getInstance()->getVFCollisionSet();
  // cout << "ees.size() = "<<ees.size() << endl;
  // cout << "vfs.size() = "<<vfs.size() << endl;

  for (size_t i = 0; i < forces.size(); ++i){
    forces[i]->vertex.clear();
    forces[i]->forces.clear();
  }

  // vertex-frace
  int vv[3] = {0};
  BOOST_FOREACH(const VFCollision &vf, vfs){

	vv[0] = _faces[vf.fid()].id0();
	vv[1] = _faces[vf.fid()].id1();
	vv[2] = _faces[vf.fid()].id2();
	const SELF_CCD::vec3f &fv0 =_verts[vv[0]];
	const SELF_CCD::vec3f &fv1 =_verts[vv[1]];
	const SELF_CCD::vec3f &fv2 =_verts[vv[2]];
	SELF_CCD::vec3f n = (fv1-fv0).cross(fv2-fv0);
	assert_gt(n.length(),0);
	n.normalize();
	
	const int v = vf.vid();
	const double deepth  = sqrt(abs((_verts[v]-_verts[vv[0]]).dot(n)));

	Vector3d f;
	f[0] = n[0];
	f[1] = n[1];
	f[2] = n[2];
	f *= (response_scalor*deepth);

	// vretex
	int obj_v, obj_id;
	convertToObejctId(v, obj_v, obj_id);
	forces[obj_id]->vertex.push_back(obj_v);
	forces[obj_id]->forces.push_back(f);

	// faces
	f = f*(-1.0f/3);
	for (int i = 0; i < 3; ++i){
	  convertToObejctId(vv[i], obj_v, obj_id);
	  forces[obj_id]->vertex.push_back(obj_v);
	  forces[obj_id]->forces.push_back(f);
	}
  }

  // // edge-edge @bug
  // int ev[4] = {0};
  // BOOST_FOREACH(const EECollision &ee, ees){

  // 	ev[0] = ee.e1v1();
  // 	ev[1] = ee.e1v2();
  // 	ev[2] = ee.e2v1();
  // 	ev[3] = ee.e2v2();

  // 	const SELF_CCD::vec3f &e1v1 =_verts[ev[0]];
  // 	const SELF_CCD::vec3f &e1v2 =_verts[ev[1]];
  // 	const SELF_CCD::vec3f &e2v1 =_verts[ev[2]];
  // 	const SELF_CCD::vec3f &e2v2 =_verts[ev[3]];

  // 	const SELF_CCD::vec3f n = ((e1v2-e1v1) - (e2v2-e2v1))*0.5f;
  // 	Vector3d f;
  // 	f[0] = n[0];
  // 	f[1] = n[1];
  // 	f[2] = n[2];
  // 	f *= (response_scalor*0.5f);

  // 	int obj_v, obj_id;
  // 	for (int i = 0; i < 2; ++i){
  // 	  convertToObejctId(ev[i], obj_v, obj_id);
  // 	  forces[obj_id]->vertex.push_back(obj_v);
  // 	  forces[obj_id]->forces.push_back(f);
  // 	}
  // 	for (int i = 2; i < 4; ++i){
  // 	  convertToObejctId(ev[i], obj_v, obj_id);
  // 	  forces[obj_id]->vertex.push_back(obj_v);
  // 	  forces[obj_id]->forces.push_back(-f);
  // 	}
  // }
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

  _EECollisionSet.clear();
  _VFCollisionSet.clear();
  ccdChecking(true);
}

void CCD_DPF::convertToObejctId(const int v, int &obj_v, int &obj_id)const{

  assert_gt(_start_index.size(),0);
  if (_start_index.size() <= 1){
	obj_v = v;
	obj_id = 0;
	return ;
  }

  if (v >= _start_index[_start_index.size()-1]){
	obj_v = v-_start_index[_start_index.size()-1];
	obj_id = _start_index.size()-1;
	return ;
  }

  for (int i = 0; i < _start_index.size()-1; ++i){
	if (_start_index[i] <= v &&  v < _start_index[i+1]){
	  obj_v = v-_start_index[i];
	  obj_id = i;
	  break;
	}
  }
}
