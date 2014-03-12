#ifndef _FULLMTLRECOVER_H_
#define _FULLMTLRECOVER_H_

#include "ComputeStiffnessMat.h"

namespace ELASTIC_OPT{
  
  /**
   * @class FullMtlRecover Given Basis W and eigenvalues lambda, tetrahedron
   * mesh, compute the materials in full space. The tetrahedron mesh also
   * contains the fixed density, and desired initial G and Lambda.
   * 
   */
  class FullMtlRecover{
	
  public:
	FullMtlRecover(){}
	
		
  protected:
	
	
  private:
	MatrixXd W;
	VectorXd Lambda;
	pTetMesh tetmesh;
	double mu_G;
	double mu_L;
  };
  
  typedef boost::shared_ptr<FullMtlRecover> pFullMtlRecover;
  
}//end of namespace

#endif /*_FULLMTLRECOVER_H_*/
