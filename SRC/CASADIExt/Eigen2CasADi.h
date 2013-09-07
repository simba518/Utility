#ifndef _EIGEN2CASADI_H_
#define _EIGEN2CASADI_H_

#include <vector>
#include <string>
#include <boost/foreach.hpp>
#include <eigen3/Eigen/Dense>
#include <symbolic/sx/sx.hpp>
#include <symbolic/matrix/matrix_tools.hpp>
#include <eigen3/Eigen/Sparse>
typedef Eigen::SparseMatrix<double> SparseMatrixD;
typedef Eigen::SparseMatrix<float> SparseMatrixF;

// #include <EigenSparseHeads.h>

/**
 * functions that convert data between eigen3 and CasADi.
 * 
 */
namespace LSW_CASADI{

  class Eigen2CasADi{
  
  public:
	static void convert(const Eigen::MatrixXd &M, CasADi::SXMatrix &SM){

	  SM.makeDense(M.rows(), M.cols(), 0.0f);
	  for (int i = 0; i < M.rows(); ++i){
		for (int j = 0; j < M.cols(); ++j){
		  SM.elem(i,j) = M(i,j);
		}
	  }
	}

	static void convert(const SparseMatrixD &M, CasADi::SXMatrix &SM){

	  std::vector<int> cols;
	  std::vector<int> rows;
	  std::vector<double> values;
	  cols.reserve(M.nonZeros());
	  rows.reserve(M.nonZeros());
	  values.reserve(M.nonZeros());
	  for (int k=0; k<M.outerSize(); ++k) {
		for (SparseMatrixD::InnerIterator it(M,k); it; ++it) {
		  cols.push_back(it.col());
		  rows.push_back(it.row());
		  values.push_back(it.value());
		}
	  }
	  SM = CasADi::DMatrix::sparse( rows, cols, values, M.rows(), M.cols() );
	}

	static void convert(const CasADi::SXMatrix &SM, SparseMatrixD &M){

	  std::vector<int> rows,cols;
	  SM.sparsity().getSparsity (rows, cols);
	  const std::vector<CasADi::SX> &data = SM.data ();

	  typedef Eigen::Triplet<double> E_Triplet;
	  std::vector<E_Triplet> SM_triplet;
	  const int nz = (int)data.size();
	  SM_triplet.reserve(nz);
	  for (int i = 0; i < nz; ++i){
		SM_triplet.push_back( E_Triplet(rows[i],cols[i],data[i].getValue()) );
	  }
	  M.resize (SM.size1(),SM.size2());
	  M.reserve (nz);
	  M.setFromTriplets( SM_triplet.begin(), SM_triplet.end() );
	}

	static void convert(const std::vector<Eigen::MatrixXd>&M,std::vector<CasADi::SXMatrix> &SM){
  
	  SM.reserve(M.size());
	  BOOST_FOREACH(const Eigen::MatrixXd &m, M){
		CasADi::SXMatrix s;
		convert(m,s);
		SM.push_back(s);
	  }
	}

	static void convert(const Eigen::VectorXd &V, CasADi::SXMatrix &SV){
  
	  SV.makeDense(V.size(),1, 0.0f);
	  for (int i = 0; i < V.size(); ++i){
		SV.elem(i) = V[i];
	  }
	}

	static void convert(const std::vector<Eigen::VectorXd>&V,std::vector<CasADi::SXMatrix> &SV){
  
	  SV.reserve(V.size());
	  BOOST_FOREACH(const Eigen::VectorXd &v, V){
		CasADi::SXMatrix s;
		convert(v,s);
		SV.push_back(s);
	  }
	}

	static CasADi::SXMatrix convert(const Eigen::VectorXd &V){

	  CasADi::SXMatrix SV;
	  SV.makeDense(V.size(),1,0.0f);
	  for (int i = 0; i < V.size(); ++i){
		SV.elem(i) = V[i];
	  }
	  return SV;
	}

	static CasADi::SXMatrix convert(const Eigen::MatrixXd &M){

	  CasADi::SXMatrix SM;
	  convert(M,SM);
	  return SM;
	}

	static Eigen::VectorXd convert2Vec(const CasADi::SXMatrix &SV){

	  Eigen::VectorXd V(SV.size());
	  for (int i = 0; i < SV.size(); ++i){
		V[i] = SV.elem(i).getValue();
	  }
	  return V;
	}

	static Eigen::MatrixXd convert(const CasADi::SXMatrix &SM){

	  Eigen::MatrixXd M(SM.size1(),SM.size2());
	  for (int i = 0; i < SM.size1(); ++i){
		for (int j = 0; j < SM.size2(); ++j){
		  M(i,j) = SM.elem(i,j).getValue();
		}
	  }
	  return M;
	}

  };

}

#endif /* _EIGEN2CASADI_H_ */
