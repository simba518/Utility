#ifndef _SPARSEGENEIGENSOLVER_H_
#define _SPARSEGENEIGENSOLVER_H_

#include <BaseSparseMatrix.h>
#include <eigensym.h>
#include <vector>
using namespace std;

namespace EIGEN3EXT{
  
  /**
   * @class SparseGenEigenSolver use aprack to solve the geneneral eigen value
   * problem K*w = lambda*M, where K and M are semi-positive definite matrix,
   * which are stored in column major sparse form. In modal analysis for
   * computing the linear modals, K is the stiffness matrix, while M is the mass
   * matrix.
   * 
   */	
  template <typename real>
  class SparseGenEigenSolver{

  public:
	static bool solve (const BaseSparseMatrix<real> &K,
					   const BaseSparseMatrix<real> &M,
					   BaseGeneralMatrix<real> &eig_vec,
					   vector<real> &eig_val,const int max_eig_num){

	  //check the parameters
	  assert(max_eig_num >0 && max_eig_num <= K.rows());
	  assert(K.rows()==M.rows()&&K.cols()==M.cols()&&K.rows()==K.cols());

	  //allocate memory
	  eig_vec.resize(K.rows(),max_eig_num);
	  eig_val.resize(max_eig_num);

	  //caculate 
	  const int k_nonzero_num = K.nonzeroNum();
	  const int k_rows = K.rows();
	  const double *k_data = &(K.Data()[0]);
	  const int *k_rowind = &(K.Rowind()[0]);
	  const int *k_colptr = &(K.Colptr()[0]);

	  const int m_nonzero_num = M.nonzeroNum();
	  const double *m_data = &(M.Data()[0]);
	  const int *m_rowind = &(M.Rowind()[0]);
	  const int *m_colptr = &(M.Colptr()[0]);

	  double *p_eig_val = &(eig_val[0]);
	  double *p_eig_vec = &(eig_vec.Data()[0]);

	  bool succ = eigensym(k_data,k_nonzero_num,k_rowind,k_colptr,
	  					   m_data,m_nonzero_num,m_rowind,m_colptr,
	  					   k_rows,max_eig_num,
						   p_eig_val,p_eig_vec);

	  return  succ;
	}
	
  };

  typedef SparseGenEigenSolver<double> SparseGenEigenSolverD;
  typedef SparseGenEigenSolver<float> SparseGenEigenSolverF;
  
}//end of namespace

#endif /*_SPARSEGENEIGENSOLVER_H_*/
