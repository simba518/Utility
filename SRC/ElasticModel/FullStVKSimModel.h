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
	FullStVKSimModel(pTetMesh tet):tetMesh(tet){}
	bool init(const std::string init_filename){
	  return false;
	}
	bool setGravity(const bool addGravity){
	  /// @todo
	  return false;
	}
	bool evaluateF(const Eigen::VectorXd &u, Eigen::VectorXd &f){
	  if(fullStvk){
		assert_eq(u.size(),rest_x.size());
		x = rest_x + u;
		fullStvk->force(x,f);
		// f *= -1.0 ???
		return true;
	  }else{
		return false;
	  }
	}
	bool evaluateK(const Eigen::VectorXd &u, SparseMatrix<double> &K_full){
	  if (fullStvk){
		assert_eq(u.size(),rest_x.size());
		x = rest_x + u;
		K_full = fullStvk->K(x);
		return true;
	  }else{
		return false;
	  }
	}
	bool evaluateK_triplet(const Eigen::VectorXd &u, VecT &K_full_t){
	  static SparseMatrix<double> K;
	  const bool succ = evaluateK(u,K);
	  getTriplet(K,K_full_t);
	  return succ;
	}
	bool evaluateM(SparseMatrix<double> &M_full){
	  if (tetMesh){
		mass.compute(M_full,*tetMesh);
		return true;
	  }else{
		return false;
	  }
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
