#include <iostream>
using namespace std;

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/CholmodSupport>
#include <eigen3/Eigen/UmfPackSupport>

#include "SparseLinearEqSolver.h"

using namespace EIGEN3EXT;

bool SparseLinearEqSolver::solve(const SparseMatrixD &K,const MatrixXd &B,MatrixXd &X,bool positive_definite){

  if(positive_definite){
	return LLTsolve(K,B,X);
  }else{
	return LUsolve(K,B,X);
  }
}

bool SparseLinearEqSolver::LLTsolve(const SparseMatrixD &K,const MatrixXd &B,MatrixXd &X){

  assert (K.rows() == K.cols());
  assert (B.rows() == K.cols());

  bool succ = true;
  X.resize(B.rows(),B.cols());
  CholmodSupernodalLLT<SparseMatrixD> solver;
  solver.compute(K);
  if(solver.info()!=Eigen::Success){
  	cout << endl<<"error: failed to decompose matrix K in solving K*X=B in SparseLinearEqSolver::LLTsolve()" << endl;
  	succ = false; // decomposiiton failed	
  }
  if(succ){
  	X = solver.solve(B);
  	if (solver.info()!=Eigen::Success){
  	  cout <<endl;
  	  cout << "error: failed to solve K*X=B in SparseLinearEqSolver::LLTsolve()" << endl;
  	  succ = false; // decomposiiton failed
  	}
  }
  return succ;
}

bool SparseLinearEqSolver::LLTsolve(const SparseMatrixD &K,const VectorXd &b,VectorXd &x){
  
  bool succ = true;
  x.resize(b.size());
  CholmodSupernodalLLT<SparseMatrixD> solver;
  solver.compute(K);
  if(solver.info()!=Eigen::Success){
  	cout << endl<<"error: failed to decompose matrix K in solving K*x=b in SparseLinearEqSolver::LLTsolve()" << endl;
  	succ = false; // decomposiiton failed	
  }
  if(succ){
  	x = solver.solve(b);
  	if (solver.info()!=Eigen::Success){
  	  cout <<endl;
  	  cout << "error: failed to solve K*x=b in SparseLinearEqSolver::LLTsolve()" << endl;
  	  succ = false; // decomposiiton failed
  	}
  }
  return succ;
}

bool SparseLinearEqSolver::LUsolve(const SparseMatrixD &K,const MatrixXd &B,MatrixXd &X){

  assert (K.rows() == K.cols());
  assert (B.rows() == K.cols());

  bool succ = true;
  X.resize(B.rows(),B.cols());
  UmfPackLU<SparseMatrixD> solver;
  solver.compute(K);
  if(solver.info()!=Eigen::Success){
  	cout << endl<<"error: failed to decompose matrix K in solving K*X=B in SparseLinearEqSolver::LUsolve()" << endl;
  	succ = false; // decomposiiton failed	
  }
  if(succ){
	X = solver.solve(B);
	if (solver.info()!=Eigen::Success){
	  cout <<endl;
	  cout << "error: failed to solve K*X=B in SparseLinearEqSolver::LUsolve()" << endl;
	  succ = false; // decomposiiton failed
	}
  }
  return succ;
}

bool SparseLinearEqSolver::LUsolve(const SparseMatrixD &K,const VectorXd &b,VectorXd &x){
  
  bool succ = true;
  UmfPackLU<SparseMatrixD> solver;
  solver.compute(K);
  if(solver.info()!=Eigen::Success){
	cout <<"failed to decompose for linear equation." << endl;
	succ = false;
  }
  if (succ) {
	x = solver.solve(b);
  }
  if(solver.info()!=Eigen::Success){
	cout << "failed to solve linear equation." <<endl;
	succ = false;
  }
  return succ;
}
