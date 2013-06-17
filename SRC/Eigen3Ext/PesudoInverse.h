#ifndef _PESUDOINVERSE_H_
#define _PESUDOINVERSE_H_

#include <eigen3/Eigen/Dense>
#include <matrixLAPACK.h>

namespace EIGEN3EXT{

#define DENSE_EIGEN_MAT Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
  
  /**
   * @class PesudoInverse compute the pesudo inverse of an eigen dense matrix.
   * 
   */
  class PesudoInverse{
	
  public:
	
	template <class T>
	static bool compute(const DENSE_EIGEN_MAT &m, DENSE_EIGEN_MAT &pinv_m,
						const T svd_threshold = 1e-6, int * rank = NULL){

	  bool succ = true;
	  pinv_m.resize(m.cols(),m.rows());
	  if(m.size() <=0 ){
		return succ;
	  }

	  const int r = m.rows();
	  const int c = m.cols();
	  T* in = const_cast<T*>(&m(0,0));
	  T* out = &pinv_m(0,0);
	  try{
		PseudoInverseMatrix(r,c,in,svd_threshold,rank,out);
	  }catch(int e){
		succ = false;
	  }
	  return succ;
	}
	
  };
  
}//end of namespace

#endif /*_PESUDOINVERSE_H_*/
