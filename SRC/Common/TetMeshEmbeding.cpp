#include <assertext.h>
#include <Log.h>
#include "TetMeshEmbeding.h"
using namespace UTILITY;

void TetMeshEmbeding::embedMesh(const VVec3d &verts){

  double w;
  VertexHash hash;

  //generate a tetrahedra mesh index
  generateHash(hash,w);
  _mapping.resize(verts.size()); ///@bug memory error in simba's note book with linux 32bit.

  //for each surface vertex of emesh find the closest tetrahedra by
  //incremental search surrounding cells.
  for (size_t i=0;i<(size_t)_mapping.size();i++) {

	//find closest tet
	const Vector3d& vert=verts[i];
	const Vector3d v(vert.x(),vert.y(),vert.z());
	size_t& tid=_mapping[i]._tid;
	Vector3d& cp=_mapping[i]._close;
	Vector4d& bc=_mapping[i]._bcCoord;
	size_t rad=0;
	while(true)	{

	  find(cp,tid,bc,rad++,v,hash,w);
	  if(tid >= 0){
		break;
	  }else{
		WARN_LOG("Cannot Find Enclosing Tet!");
	  }
	}
  }
}

// void TetMeshEmbeding::find(Vector3d& cp,size_t& tid,Vector4d& bary,const size_t& rad,
// 						  const Vector3d& pos,const VertexHash& hash,const double& w)const{

//   boost::unordered_set<size_t> tested;
//   VertexHash::const_iterator iter;
//   double thres=w*(rad+0.5f);thres*=thres;
//   double sqrMinDist=-1.0f,sqrDist;
//   tid=-1;
//   Vector3d cpTmp(0,0,0);
//   Vector4d baryTmp;
//   const Vector3i cid=floor((Vector3i)(pos/w));
//   for(size_t x=cid.x()-rad;x<=cid.x()+rad;x++)
// 	for(size_t y=cid.y()-rad;y<=cid.y()+rad;y++)
// 	  for(size_t z=cid.z()-rad;z<=cid.z()+rad;z++)
// 		{
// 		  if((iter=hash.find(HashedId(x,y,z,0))) != hash.end())
// 			{
// 			  for(std::vector<size_t>::const_iterator beg=iter->second.begin(),end=iter->second.end();beg!=end;beg++)
// 				{
// 				  if(tested.find(*beg) == tested.end())
// 					{
// 					  const Vector4i& tet=_vol_mesh->tets()[*beg];
// 					  tetrahedron t(_vol_mesh->nodes()[tet[0]],_vol_mesh->nodes()[tet[1]],
// 									_vol_mesh->nodes()[tet[2]],_vol_mesh->nodes()[tet[3]]);
// 					  t.calcPointDist(pos,sqrDist,cpTmp,baryTmp);
// 					  if(sqrDist < thres && (sqrMinDist < 0.0f || sqrDist < sqrMinDist))
// 						{
// 						  sqrMinDist=sqrDist;
// 						  tid=*beg;
// 						  bary=t.bary(pos);
// 						  cp=cpTmp;
// 						  if(sqrMinDist == 0.0f)
// 							return;
// 						}
// 					  tested.insert(*beg);
// 					}
// 				}
// 			}
// 		}
// }

// void TetMeshEmbeding::generateHash(VertexHash& hash,double& w){

//   //find the tet is minimal volume
//   double vol=_vol_mesh->volume(0);
//   for(size_t i=1;i<(size_t)_vol_mesh->tets().size();i++)
// 	vol+=_vol_mesh->volume(i);
//   //approximate hash width
//   w=std::pow(vol/(double)_vol_mesh->tets().size(),(double)1.0f/3.0f);
//   //insert all tets
//   hash.clear();
//   for(size_t i=0;i<(size_t)_vol_mesh->tets().size();i++)
// 	{
// 	  //calculate bounding box
// 	  const Vector4i& tet=_vol_mesh->tets()[i];
// 	  bBox<double> bb;
// 	  bb.setPoints(_vol_mesh->nodes()[tet[0]],_vol_mesh->nodes()[tet[1]],_vol_mesh->nodes()[tet[2]]);
// 	  bb.setUnion(_vol_mesh->nodes()[tet[3]]);
// 	  //insert to all cells
// 	  const Vector3di minC=floor((Vector3d)(bb._minC/w));
// 	  const Vector3di maxC=ceil((Vector3d)(bb._maxC/w));
// 	  for(size_t x=minC.x();x<=maxC.x();x++)
// 		for(size_t y=minC.y();y<=maxC.y();y++)
// 		  for(size_t z=minC.z();z<=maxC.z();z++)
// 			hash[HashedId(x,y,z,0)].push_back(i);
// 	}
// }

void TetMeshEmbeding::generateNewVertexList(const VectorXd& nodes,VVec3d& verts){

  verts.resize(_mapping.size());
  for(size_t i=0;i<(size_t)_mapping.size();i++) {

	  const VertexMap& vm=_mapping[i];
	  const Vector4i& tet=_vol_mesh->tets()[vm._tid];
	  Vector3d& pt=verts[i];
	  pt= nodes.block<3,1>(tet[0]*3,0)*vm._bcCoord.x()+
		nodes.block<3,1>(tet[1]*3,0)*vm._bcCoord.y()+
		nodes.block<3,1>(tet[2]*3,0)*vm._bcCoord.z()+
		nodes.block<3,1>(tet[3]*3,0)*vm._bcCoord.w();
	}
}

void TetMeshEmbeding::interpolateForces(const vector<int>& node_ids,const VVec3d& forces, VectorXd& rlst_forces){

  /// @todo write testcases for void TetMeshEmbeding::interpolateForces
  assert_eq (rlst_forces.size()%3,0);
  assert_eq (node_ids.size(),forces.size());
  for (size_t i = 0; i < node_ids.size(); ++i){

	const VertexMap& vm=_mapping[node_ids[i]];
	const Vector4i& tet=_vol_mesh->tets()[vm._tid];

	rlst_forces.block<3,1>(tet[0]*3,0) += vm._bcCoord.x()*forces[i];
	rlst_forces.block<3,1>(tet[1]*3,0) += vm._bcCoord.y()*forces[i];
	rlst_forces.block<3,1>(tet[2]*3,0) += vm._bcCoord.z()*forces[i];
	rlst_forces.block<3,1>(tet[3]*3,0) += vm._bcCoord.w()*forces[i];
  }
}
