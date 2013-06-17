#ifndef _EIGENSPAREMASTRIXTOOLS_H_
#define _EIGENSPAREMASTRIXTOOLS_H_

#include <vector>
#include <set>
#include <eigen3/Eigen/Dense>
#include <EigenSparseHeads.h>
#include <assertext.h>
#include <boost/foreach.hpp>

namespace EIGEN3EXT{

  /************************************ create ********************************/
  template <class T>  
  const Eigen::SparseMatrix<T> &createFromDense(const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> &M, 
												Eigen::SparseMatrix<T> &S, const T tol = 1e-16){

	typedef Eigen::Triplet<T> E_Triplet;
	std::vector<E_Triplet> striplet;
	striplet.reserve(M.size());
	for (int i = 0; i < M.rows(); ++i) {
	  for (int j = 0; j < M.cols(); ++j) {
		if ( fabs(M(i,j)) >= tol )
		  striplet.push_back( E_Triplet(i,j,M(i,j)) );
	  }
	}
	S.resize(M.rows(), M.cols());
	S.setFromTriplets(striplet.begin(), striplet.end());
	return S;
  }

  template <class T> 
  const Eigen::SparseMatrix<T> createFromDense(const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> &M,
											   const T tol = 1e-16){
	
	Eigen::SparseMatrix<T> S;
	return createFromDense(M,S,tol);
  }

  template <class T>
  const Eigen::SparseMatrix<T> &eye(const int n, Eigen::SparseMatrix<T> &A, const T value){

	assert_ge(n,0);
	typedef Eigen::Triplet<T> Tri;
	std::vector<Tri> triplets;
	triplets.reserve(n);
	for (int i = 0; i < n; ++i){
	  triplets.push_back( Tri(i, i, value) );
	}
	A.resize(n,n);
	A.reserve(n);
	A.setFromTriplets( triplets.begin(),triplets.end() );
	A.makeCompressed();
	return A;
  }

  template <class T>
  const Eigen::SparseMatrix<T> eye(const int n, const T value){

	assert_ge(n,0);
	Eigen::SparseMatrix<T> S;
	return eye(n,S,value);
  }

  template <class T, class MatrixType> 
  const Eigen::SparseMatrix<T> &eye(Eigen::SparseMatrix<T> &A_block_diag, const std::vector<MatrixType> &block_mats){

	const size_t mat_num = block_mats.size();
	if (mat_num <=0 ){
	  A_block_diag.resize(0,0);
	  return A_block_diag;
	}

	const size_t n = block_mats[0].size()*mat_num;
	typedef Eigen::Triplet<T> Tri;
	std::vector<Tri> triplets;
	triplets.reserve(n);

	int rows = 0;
	int cols = 0;
	for (size_t mi = 0; mi < mat_num; ++mi){
	  const int r = block_mats[mi].rows();
	  const int c = block_mats[mi].cols();
	  for (int i = 0; i < r; ++i){
		for (int j = 0; j < c; ++j){
		  triplets.push_back( Tri(i+rows, j+cols, block_mats[mi](i,j)) );
		}
	  }
	  rows += r;
	  cols += c;
	}

	A_block_diag.resize(rows,cols);
	A_block_diag.reserve(triplets.size());
	A_block_diag.setFromTriplets( triplets.begin(),triplets.end() );
	A_block_diag.makeCompressed();
	return A_block_diag;
  }

  template <class T> 
  const Eigen::SparseMatrix<T> random(const int r, const int c, const T scalor){

	Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic> M = Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic>::Random(r,c);
	M *= scalor;
	Eigen::SparseMatrix<T> S;
	return createFromDense(M,S);
  }


  
  /************************************ block *********************************/

  /** 
   * Generate the sparse matrix P which will remove the rows of Matrix A in
   * remove_rows_set by using:
   * A' = P*A.
   * 
   * @param total_rows rows of the original matrix A.
   * @param remove_rows_set the rows to be remove.
   * @param P output matrix.
   * @param remove if it is false, the rows in remove_rows_set will be
   * preserved, and others will be removed.
   * @note all indices in remove_rows_set should be in [0,total_rows-1].
   */	
  template <class T>
  const Eigen::SparseMatrix<T> &genReshapeMatrix(const int total_rows, 
												 const std::set<int> &remove_rows_set, 
												 Eigen::SparseMatrix<T> &P,
												 const bool remove = true){
	  
	typedef Eigen::Triplet<T> Tri;

	const int rows = remove ? total_rows - (int)remove_rows_set.size():(int)remove_rows_set.size();
	const int cols = total_rows;
	const int nonzeros = rows;

	std::vector<Tri> P_triplets;
	P_triplets.reserve(nonzeros);
	  
	if(remove){
	  for (int i = 0; i < total_rows; ++i){
		if ( remove_rows_set.find(i) == remove_rows_set.end() ){
		  P_triplets.push_back( Tri((int)P_triplets.size(), i, 1) );
		}
	  }
	}else{
	  for (int i = 0; i < total_rows; ++i){
		if ( remove_rows_set.find(i) != remove_rows_set.end() ){
		  P_triplets.push_back( Tri((int)P_triplets.size(), i, 1) );
		}
	  }
	}

	P.resize(rows,cols);
	P.reserve( nonzeros );
	P.setFromTriplets( P_triplets.begin(),P_triplets.end() );
	P.makeCompressed();
	return P;
  }

  template <class T>
  Eigen::SparseMatrix<T> genReshapeMatrix(const int total_rows, 
										  const std::set<int> &remove_rows_set, 							 
										  const bool remove = true){
	Eigen::SparseMatrix<T> P;
	genReshapeMatrix(total_rows,remove,P,remove);
	return P;
  }

  /** 
   * Generate a P that will remove i-th r sub-rows from Matrix A,
   * A' = P*A.
   * example:
   *     |r1|                                          
   * A = |r2|, r = 2, remove_rows_set = {0}, then A' = |r3|.
   *     |r3|
   */
  template <class T>
  const Eigen::SparseMatrix<T> &genReshapeMatrix(const int total_rows, 
												 const int r,
												 const std::set<int> &remove_rows_set,
												 Eigen::SparseMatrix<T> &P, 
												 const bool remove = true){
	  
	std::set<int> rm_rows_set;
	BOOST_FOREACH(int ele, remove_rows_set){
	  if (ele*r >= 0 && ele*r + r <= total_rows){
		for (int i = 0; i < r; ++i){
		  rm_rows_set.insert(ele*r + i);
		}
	  }
	}
	genReshapeMatrix(total_rows,rm_rows_set,P, remove);
	return P;
  }

  template <class T>
  Eigen::SparseMatrix<T> genReshapeMatrix(const int total_rows, 
										  const int r,
										  const std::set<int> &remove_rows_set,							    
										  const bool remove = true){
	Eigen::SparseMatrix<T> P;
	genReshapeMatrix(total_rows,r,remove_rows_set,P, remove);
	return P;
  }

  template <class T> 
  const Eigen::SparseMatrix<T> block(const Eigen::SparseMatrix<T> &S, 
									 const int r0, const int c0, 
									 const int rows, const int cols){
	
	assert_ge (r0,0);
	assert_ge (c0,0);
	assert_gt (rows,0);
	assert_gt (cols,0);
	assert_le (r0+rows, S.rows());
	assert_le (c0+cols, S.cols());

	std::set<int> remove_rows_set, remove_cols_set;
	for (int i = 0; i < rows; ++i)  remove_rows_set.insert(r0+i);
	for (int i = 0; i < cols; ++i)  remove_cols_set.insert(c0+i);

	Eigen::SparseMatrix<T> B, P1, P2;
	genReshapeMatrix(S.rows(), remove_rows_set, P1, false);
	genReshapeMatrix(S.cols(), remove_cols_set, P2, false);
	assert_eq (P1.cols(), S.rows());
	assert_eq (P2.cols(), S.cols());

	B = (P1*S)*(P2.transpose());
	return B;
  }


  /************************************ inverse *********************************/
  template <class T>
  const Eigen::SparseMatrix<T> &inverse(const Eigen::SparseMatrix<T> &P, Eigen::SparseMatrix<T> &inv_P){

	/// @todo function inverse(), very slow.
	Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic> M = P;
	M = M.inverse();
	createFromDense(M, inv_P);
	return inv_P;
  }

  template <class T>
  const Eigen::SparseMatrix<T> inverse(const Eigen::SparseMatrix<T> &P){

	Eigen::SparseMatrix<T> inv_P;
	return inverse(P,inv_P);
  }

  

}
#endif /* _EIGENSPAREMASTRIXTOOLS_H_ */
