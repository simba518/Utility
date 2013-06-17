#ifndef _SPARSELINEAREQSOLVER_H_
#define _SPARSELINEAREQSOLVER_H_

#include <eigen3/Eigen/Dense>
using namespace Eigen;

#include <EigenSparseHeads.h>

namespace EIGEN3EXT{
  
  /**
   * @class SparseLinearEqSolver solve the linear siystem K*X = B where K is a
   * sparse matrix, X and B are dense matrix (vectors).
   * @see http://eigen.tuxfamily.org/dox-devel/TutorialSparse.html
   */
  class SparseLinearEqSolver{
	
  public:
	// solve K*B = X, using llt (if K is positive_definite) or lu.
	static bool solve(const SparseMatrixD &K,const MatrixXd &B,MatrixXd &X,bool positive_definite);

	// solve K*B = X using llt ,where K should be positive_definite.
	static bool LLTsolve(const SparseMatrixD &K,const MatrixXd &B,MatrixXd &X);

	// solve K*b = x using llt ,where K should be positive_definite.
	static bool LLTsolve(const SparseMatrixD &K,const VectorXd &b,VectorXd &x);

	// solve K*B = X using lu
	static bool LUsolve(const SparseMatrixD &K,const MatrixXd &B,MatrixXd &X);

	// solve K*b = x using lu
	static bool LUsolve(const SparseMatrixD &K,const VectorXd &b,VectorXd &x);
	
  };
  
}//end of namespace

#endif /*_SPARSELINEAREQSOLVER_H_*/
