#ifndef _FULLELASTICMODEL_H_
#define _FULLELASTICMODEL_H_

#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <ElasticForceTetFullStVK.h>
#include <TetMesh.h>
#include <MassMatrix.h>
using namespace std;
using namespace Eigen;
using namespace UTILITY;

namespace SIMULATOR{

  typedef vector<Eigen::Triplet<double> > VecT;

  /**
   * @class BaseFullModel interface for computing internal forces, stiffness
   * matrix, mass matrix in full space, which are used in full simulator.
   */
  class BaseFullModel{
	
  public:
	// initialize from file.
	virtual bool init(const std::string init_filename) = 0;
	virtual bool prepare() = 0;

	// compute the internal forces.
	virtual bool evaluateF(const Eigen::VectorXd &u, Eigen::VectorXd &f) = 0;

	// compute the stiffness matrix or its triplet. both in full format.
	virtual bool evaluateK(const Eigen::VectorXd &u, SparseMatrix<double> &K_full)=0;
	virtual bool evaluateK_triplet(const Eigen::VectorXd &u, VecT &K_full_t)=0;

	// compute the mass matrix or its triplet. both in full format.
	virtual bool evaluateM(SparseMatrix<double> &M_full)=0;
	virtual bool evaluateM_triplet(VecT &M_full_t)=0;

	// dimension of the full space.
	virtual int dimension()const = 0;
	virtual const VectorXd &getRestShape()const = 0;
  };
  typedef boost::shared_ptr<BaseFullModel> pBaseFullModel;

  /**
   * @class FullStVKSimModel data model for full stvk simulation.
   * 
   */
  class FullStVKSimModel:public BaseFullModel{
	
  public:
	FullStVKSimModel(){
	  fullStvk = pElasticForceTetFullStVK (new ElasticForceTetFullStVK());
	}
	bool init(const std::string filename){return true;}
	bool prepare(){
	  
	  if (tetMesh){
		tetMesh->nodes(rest_x);
		x = rest_x;
		fullStvk->prepare();
		return true;
	  }
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
	const VectorXd &getRestShape()const{
	  return rest_x;
	}

	void setTetMesh(pTetMesh_const tet){
	  tetMesh = tet;
	  fullStvk->setVolMesh(tetMesh);
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
	pTetMesh_const tetMesh;
	pElasticForceTetFullStVK fullStvk;
	MassMatrix mass;
	VectorXd rest_x;
	VectorXd x;
  };
  typedef boost::shared_ptr<FullStVKSimModel> pFullStVKSimModel;
  
}//end of namespace

#endif /* _FULLELASTICMODEL_H_ */
