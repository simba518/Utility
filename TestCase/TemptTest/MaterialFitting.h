#ifndef _MATERIALFITTING_H_
#define _MATERIALFITTING_H_

#include <iostream>
#include <MatrixTools.h>
#include <MatrixIO.h>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <CASADITools.h>
#include <MassMatrix.h>
#include <ElasticForceTetFullStVK.h>
#include <Timer.h>
#include "ComputeStiffnessMat.h"
#include "ComputeMassMat.h"
using namespace std;
using namespace UTILITY;
using namespace EIGEN3EXT;
using namespace ELASTIC_OPT;
using namespace CASADI;
using CasADi::SX;

namespace ELASTIC_OPT{
  
  /**
   * @class MaterialFitting
   * 
   */
  class MaterialFitting{
	
  public:
	MaterialFitting(){
	  tetmesh = pTetMesh(new TetMesh());
	  use_hessian = false;
	  mu_neigh_G = 1.0f;
	  mu_neigh_L = 1.0f;
	  mu_neigh_rho = 1.0f;
	  mu_average_rho = 1.0f;
	}
	void loadTetMesh(const string filename){
	  const bool succ = tetmesh->load(filename); assert(succ);
	}
	void loadMtl(const string filename){
	  bool succ = tetmesh->loadElasticMtl(filename); assert(succ);
	}
	void loadFixednodes(const string filename){
	  vector<int> fixednodes_vec;
	  const bool succ = loadVec(filename,fixednodes_vec, TEXT); assert(succ);
	  fixednodes.clear();
	  for (int i = 0; i < fixednodes_vec.size(); ++i){
		fixednodes.insert(fixednodes_vec[i]);
	  }
	}
	void setWLambda(const MatrixXd &eigen_W,const VectorXd &eigen_lambda){
	  CASADI::convert(eigen_W, W);
	  lambda.resize(eigen_lambda.size(), eigen_lambda.size());
	  for (int i = 0; i < eigen_lambda.size(); ++i)
		lambda.elem(i,i) = eigen_lambda[i];
	  assert_lt(W.size2(), W.size1());
	  assert_eq(lambda.size1(), W.size2());
	  assert_eq(lambda.size2(), W.size2());
	}
	void useHessian(const bool use){
	  use_hessian = use;
	}
	void setMuSmooth(const double mu_G,const double mu_L,const double mu_rho){
	  assert_ge(mu_G,0.0f);
	  assert_ge(mu_L,0.0f);
	  assert_ge(mu_rho,0.0f);
	  mu_neigh_G = mu_G;
	  mu_neigh_L = mu_L;
	  mu_neigh_rho = mu_rho;
	}
	void setMuAverageDensity(const double mu){
	  assert_ge(mu,0.0f);
	  mu_average_rho = mu;
	}
	void computeK();
	void computeM();
	void removeFixedDOFs();
	virtual void assembleObjfun();
	virtual void solve();
	virtual void saveResults(const string filename)const;
	
  protected:
	virtual void initShearG(const int num_tet, VSX &G)const{
	  G = CASADI::makeSymbolic(num_tet, "G");
	}
	virtual void initLame(const int num_tet, VSX &Lame)const{
	  Lame = CASADI::makeSymbolic(num_tet, "lame");
	}
	virtual void initDensity(const int num_tet, VSX &rho)const{
	  rho = CASADI::makeSymbolic(num_tet, "rho");
	}
	virtual void initAllVariables(VSX &x)const{
	  x = CASADI::connect(CASADI::connect(G,Lame),rho);
	}
	virtual void getInitValue(VectorXd &init_x)const;
	virtual vector<double> getShearGResult()const;
	virtual vector<double> getLameResult()const;
	virtual vector<double> getDensityResult()const;
	void addSmoothObjfun(SX &objfun)const;
	void addAverageDensityObjfun(SX &objfun)const;
	void addFixedNodesObjfun(SX &objfun)const;
	
  protected:
	SXMatrix W;
	SXMatrix lambda;
	pTetMesh tetmesh;
	set<int> fixednodes;
	bool use_hessian;
	double mu_neigh_G, mu_neigh_L, mu_neigh_rho, mu_average_rho;
	SX objfun;
	CASADI::VSX G;
	CASADI::VSX Lame;
	CASADI::VSX rho;
	SXMatrix K;
	SXMatrix M;
	vector<double> rlst;
  };
  typedef boost::shared_ptr<MaterialFitting> pMaterialFitting;
  
  class MaterialFittingM2: public MaterialFitting{

  public:
	MaterialFittingM2():MaterialFitting(){
	  mu_mass = 1.0f;
	  mu_stiff = 1.0f;
	}
	void setMuMass(const double mu){
	  assert_ge(mu,0.0f);
	  mu_mass = mu;
	}
	void setMuStiff(const double mu){
	  assert_ge(mu,0.0f);
	  mu_mass = mu;
	}
	void assembleObjfun();

  private:
	double mu_mass;
	double mu_stiff;
  };

  class MaterialFittingM3: public MaterialFitting{

  public:
	MaterialFittingM3():MaterialFitting(){
	  mu_mass = 1.0f;
	  mu_stiff = 1.0f;
	}

  protected:
	void initDensity(const int num_tet, VSX &rho)const;
	void initAllVariables(VSX &x)const{
	  x = CASADI::connect(G,Lame);
	}
	void getInitValue(VectorXd &init_x)const;
	vector<double> getDensityResult()const{
	  return tetmesh->material()._rho;
	}

  private:
	double mu_mass;
	double mu_stiff;
  };

}//end of namespace

#endif /*_MATERIALFITTING_H_*/
