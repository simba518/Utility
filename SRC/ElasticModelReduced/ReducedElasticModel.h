#ifndef _REDUCEDELASTICMODEL_H_
#define _REDUCEDELASTICMODEL_H_

#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <eigen3/Eigen/Dense>
#include <assertext.h>
using namespace std;
using namespace Eigen;

namespace SIMULATOR{
  
  /**
   * @class ReducedElasticModel base class for reduced elastic data model.
   * 
   */
  class ReducedElasticModel{
	
  public:
	// initialize from the ini file
	virtual bool init(const string init_filename) = 0;

	// compute the internal forces in subspace
	virtual bool evaluateF(const VectorXd &reduced_u,VectorXd &f) = 0;

	// compute the stiffness matrix in subspace
	virtual bool evaluateK(const VectorXd &reduced_u,MatrixXd &K) = 0;

	// return the reduced mass matrix: U^t*M*U
	virtual const MatrixXd &getReducedMassM()const{
	  assert_eq(M.rows(), M.cols());
	  assert_eq(M.rows(), reducedDim());
	  return M;
	}

	// convert the displacements in subsapce to fullspace
	virtual bool computeFullDisp(const VectorXd &reduced_u,VectorXd &full_u){
	  assert_eq(reduced_u.size(), B.cols());
	  full_u = B*reduced_u;
	  return true;
	}

	// dimension of the reduced subspace
	int reducedDim()const{
	  return B.cols();
	}

	// dimension of the full space.
	int fullDim()const{
	  return B.rows();
	}

	// get basis of the subspace
	const MatrixXd &getModalBasis()const{
	  return B;
	}
	
  protected:
	MatrixXd B;
	MatrixXd M;
  };
  typedef boost::shared_ptr<ReducedElasticModel> pReducedElasticModel;

  /**
   * @class CubaturedElasticModel base class for reduced elastic data model.
   * 
   */
  class CubaturedElasticModel: public ReducedElasticModel{
	
  public:
	void setCubature(const vector<double> &w, const vector<int> &S){
	  assert_eq(w.size(), S.size());
	  weights = w;
	  sampledTets = S;
	}
	bool init(const string init_filename);
	bool evaluateF(const VectorXd &reduced_u,VectorXd &f);
	bool evaluateK(const VectorXd &reduced_u,MatrixXd &K);

  private:
	vector<double> weights;
	vector<int> sampledTets;
  };
  typedef boost::shared_ptr<CubaturedElasticModel> pCubaturedElasticModel;
  
}//end of namespace

#endif /*_REDUCEDELASTICMODEL_H_*/
