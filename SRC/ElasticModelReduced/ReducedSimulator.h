#ifndef _REDUCEDSIMULATOR_H_
#define _REDUCEDSIMULATOR_H_

#include <set>
#include <ReducedElasticModel.h>
using namespace std;

namespace SIMULATOR{
  
  /**
   * @class ReducedSimulator base class for subspace simulation.
   * 
   */
  class ReducedSimulator{
	
  public:
	// set
	ReducedSimulator(pReducedElasticModel m):model(m){
	  h = 1.0f;
	  alpha_k = 0.0f;
	  alpha_m = 0.0f;
	  reset();
	}
	virtual bool init(const string init_filename);
	void setQ0(const VectorXd &q0){q = q0;}
	void setV0(const VectorXd &v0){v = v0;}
	void setTimeStep(const double h){
	  assert_gt (h, 0.0f);
	  this->h = h;  
	}
	void setDampings(const double alpha_k,const double alpha_m){
	  assert_ge (alpha_k, 0.0f);
	  assert_ge (alpha_m, 0.0f);
	  this->alpha_k = alpha_k;
	  this->alpha_m = alpha_m;  
	}
	virtual void reset(){
	  const int r = reducedDim();
	  v.resize(r);
	  v.setZero();
	  q.resize(r);
	  q.setZero();
	  red_fext.resize(r);
	  red_fext.setZero();
	}

	// set constraints
	virtual void setConGroups(const vector<set<int> > &groups) = 0;
	virtual void setUc(const VectorXd &uc) = 0;
	virtual void removeAllCon() = 0;

	// set external forces in full space
	virtual void setExtForce(const VectorXd &full_ext){
	  assert(model);
	  const MatrixXd &B = model->getModalBasis();
	  assert_eq(B.rows(), full_ext.size());
	  red_fext = B.transpose()*full_ext;
	}

	// simulate
	virtual bool forward() = 0;
	
	// get
	const VectorXd &getQ()const{return q;}
	const VectorXd &getV()const{return v;}
	double getTimestep()const{
	  return h;
	}
	double getAlphaK()const{
	  return alpha_k;
	}
	double getAlphaM()const{
	  return alpha_m;
	}
	int reducedDim()const{
	  return NULL==model ? model->reducedDim():0;
	}
	int fullDim()const{
	  return NULL==model ? model->fullDim():0;
	}

  protected:
	pReducedElasticModel model;
	double h;
	double alpha_k;
	double alpha_m;

	VectorXd v; // velocity.
	VectorXd q; // position in sub space.
	VectorXd red_fext; // external forces in subspace.
  };
  typedef boost::shared_ptr<ReducedSimulator> pReducedSimulator;

  /**
   * @class ReducedImpLogConSimulator implicit integration with Lagrangian constraints.
   * 
   */  
  class ReducedImpLogConSimulator:public ReducedSimulator{
	
  public:
	ReducedImpLogConSimulator(pReducedElasticModel m):ReducedSimulator(m){}
	void setConGroups(const vector<set<int> > &groups);
	void setUc(const VectorXd &uc){
	  this->uc = uc;
	}
	void removeAllCon();
	bool forward();

  private:
	MatrixXd A;
	VectorXd b;
	MatrixXd C;
	VectorXd uc;
  };
  
}//end of namespace

#endif /*_REDUCEDSIMULATOR_H_*/
