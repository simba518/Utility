#ifndef _FULLSTVKSIMMODEL_H_
#define _FULLSTVKSIMMODEL_H_

#include <boost/shared_ptr.hpp>
#include <TetMesh.h>
#include <BaseFullModel.h>
#include <ElasticForceTetFullStVK.h>
#include <MassMatrix.h>
using namespace UTILITY;

namespace SIMULATOR{
  
  /**
   * @class FullStVKSimModel data model for full stvk simulation.
   * 
   */
  class FullStVKSimModel:public BaseFullModel{
	
  public:
	FullStVKSimModel(){
	  tetMesh = pTetMesh(new TetMesh);
	}
	FullStVKSimModel(pTetMesh tet):tetMesh(tet){}
	bool init(const std::string init_filename){
	  
	  if (!tetMesh){
		tetMesh = pTetMesh(new TetMesh);
	  }

	  bool succ = false;
	  JsonFilePaser jsonf;
	  if (jsonf.open(init_filename)){
		string vol_file, elastic_mtl;
		succ = jsonf.readFilePath("vol_file", vol_file);
		succ &= jsonf.readFilePath("elastic_mtl", elastic_mtl);
		succ &= tetMesh->load(vol_file);
		succ &= tetMesh->loadElasticMtl(elastic_mtl);
	  }
	  if (succ){
		fullStvk = pElasticForceTetFullStVK (new ElasticForceTetFullStVK(tetMesh));
	  }

	  tetMesh->nodes(rest_x);
	  x = rest_x;
	  ERROR_LOG_COND("failed to initialze. ", succ);
	  return succ;
	}
	bool setGravity(const bool addGravity){
	  /// @todo
	  assert(false);
	  return false;
	}
	bool evaluateF(const Eigen::VectorXd &u, Eigen::VectorXd &f){
	  assert(fullStvk);
	  assert_eq(u.size(),rest_x.size());
	  x = rest_x + u;
	  fullStvk->force(x,f);
	  f = -f; /// @todo
	  return true;
	}
	bool evaluateK(const Eigen::VectorXd &u, SparseMatrix<double> &K_full){
	  assert(fullStvk);
	  assert_eq(u.size(),rest_x.size());
	  x = rest_x + u;
	  K_full = fullStvk->K(x);
	  K_full *= -1.0f;
	  return true;
	}
	bool evaluateK_triplet(const Eigen::VectorXd &u, VecT &K_full_t){
	  static SparseMatrix<double> K;
	  const bool succ = evaluateK(u,K);
	  getTriplet(K,K_full_t);
	  return succ;
	}
	bool evaluateM(SparseMatrix<double> &M_full){
	  assert(tetMesh);
	  mass.compute(M_full,*tetMesh);
	  return true;
	}
	bool evaluateM_triplet(VecT &M_full_t){
	  static SparseMatrix<double> M;
	  const bool succ = evaluateM(M);
	  getTriplet(M,M_full_t);
	  return succ;
	}
	int dimension()const{
	  return rest_x.size();
	}
	pTetMesh_const getTetMesh()const{
	  return tetMesh;
	}

  protected:
	static void getTriplet(const SparseMatrix<double> &mat, VecT&trip){
	  trip.clear();
	  trip.reserve(mat.nonZeros());
	  for (int k=0; k <mat.outerSize(); ++k)
		for (SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
		  trip.push_back(Triplet<double>(it.row(),it.col(),it.value()));
	}
	
  private:
	pTetMesh tetMesh;
	pElasticForceTetFullStVK fullStvk;
	MassMatrix mass;
	VectorXd rest_x;
	VectorXd x;
  };
  
  typedef boost::shared_ptr<FullStVKSimModel> pFullStVKSimModel;
  
}//end of namespace

#endif /*_FULLSTVKSIMMODEL_H_*/
