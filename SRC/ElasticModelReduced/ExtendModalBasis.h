#ifndef _EXTENDMODALBASIS_H_
#define _EXTENDMODALBASIS_H_

#include <eigen3/Eigen/Dense>
using namespace Eigen;

namespace SIMULATOR{
  
  /**
   * @class ExtendModalBasis generate the extended modal basis.
   * @see An efficient construction of reduced deformable objects, siggraph asia, 2013.
   */
  class ExtendModalBasis{
	
  public:
	template<class T>
	static void construct(Matrix<T,-1,-1> &linearBasis, Matrix<T,-1,-1> &extendedBasis){
	  
	}

	template<class T>
	static void orthonormalize(Matrix<T,-1,-1> &extendedBasis){
	  
	}
	
  };
  
}//end of namespace

#endif /*_EXTENDMODALBASIS_H_*/
