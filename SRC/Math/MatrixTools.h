#ifndef _MATRIXTOOLS_H_
#define _MATRIXTOOLS_H_

#include <vector>
#include <set>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SVD>
#include <assertext.h>
#include <boost/foreach.hpp>

namespace EIGEN3EXT{

  template <class T>
  struct MatType{
    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mat;
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1> vector;
  };

  /************************************ create ********************************/
  template <class T>
  const Eigen::Matrix<T,3,3> &createFromRowMajor (Eigen::Matrix<T,3,3> &M, const T *A){

	const int r = 3;
	const int c = 3;
	for (int i = 0; i < r; ++i) {
	  for (int j = 0; j < c; ++j) {
		M(i,j) = A[i*c + j];
	  }
	}
	return M;
  }


  template <class T>
  const typename MatType<T>::mat &createFromRowMajor(typename MatType<T>::mat &M,
													 const T *A, const int r, const int c){

	assert_ge (r,0);
	assert_ge (c,0);
	M.resize(r,c);
	for (int i = 0; i < r; ++i) {
	  for (int j = 0; j < c; ++j) {
		M(i,j) = A[i*c + j];
	  }
	}
	return M;
  }

  template <class T>
  void createRowMajor(const typename MatType<T>::mat &M,T *A){

	const int r = M.rows();
	const int c = M.cols();
	for (int i = 0; i < r; ++i) {
	  for (int j = 0; j < c; ++j) {
		A[i*c + j] = M(i,j);
	  }
	}
  }

  template <class T>
  const typename MatType<T>::mat &createFromColMajor(typename MatType<T>::mat &M, 
  													 const T *A, const int r, const int c){
  	assert_ge (r,0);
  	assert_ge (c,0);
  	M.resize(r,c);
  	for (int i = 0; i < r; ++i) {
  	  for (int j = 0; j < c; ++j) {
  		M(i,j) = A[j*r + i];
  	  }
  	}
  	return M;
  }

  /************************************ decomposition *************************/
  
  /*
   * Modified Sigular Value Decomposition.
   * F = U*D*Vt, where U and Vt are pure rotation matrices, e.g. det(U)=det(V)=1.
   * @see study record: Polar Decomposition.
   */
  template <class T>
  void ModifiedSVD3x3(const Eigen::Matrix<T,3,3> &F,
					  Eigen::Matrix<T,3,3> &U,
					  Eigen::Matrix<T,3,3> &Vt,
					  Eigen::Matrix<T,3,3> &D){

	Eigen::JacobiSVD<Eigen::Matrix<T,3,3> > svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
	U = svd.matrixU();
	Vt = svd.matrixV().transpose();
	assert_eq (svd.singularValues().size(),3);
	D.setZero();
	D(0,0) = svd.singularValues()[0];
	D(1,1) = svd.singularValues()[1];
	D(2,2) = svd.singularValues()[2];

	if (U.determinant() < 0) {
	  U.col(2) *= (-1.0f);
	  D(2,2) *= -1.0f;
	}
	if (Vt.determinant() < 0){
	  Vt.row(2) *= (-1.0f);
	  D(2,2) *= -1.0f;
	}
  }

  /*
   * Modified Polar Decomposition.
   * F = RS, where R is a pure rotation matrix, e.g. det(R)=1, and S is symetric
   * but maybe not positive-definit.
   * @see study record: Polar Decomposition.
   */
  template <class T>
  void ModifiedPD3x3(const Eigen::Matrix<T,3,3> &F,
  					 Eigen::Matrix<T,3,3> &R,
  					 Eigen::Matrix<T,3,3> &S){

  	Eigen::Matrix<T,3,3> &U = R;
	Eigen::Matrix<T,3,3> &Vt = S;
	Eigen::Matrix<T,3,3> D;
  	ModifiedSVD3x3(F,U,Vt,D);
  	R = U*Vt;
  	S = Vt.transpose()*D*Vt;
  }
  
}

#endif /* _MATRIXTOOLS_H_ */
