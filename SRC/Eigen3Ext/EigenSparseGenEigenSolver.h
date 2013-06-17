#ifndef _EIGENSPARSEGENEIGENSOLVER_H_
#define _EIGENSPARSEGENEIGENSOLVER_H_

#include <eigen3/Eigen/Dense>
#include <EigenSparseHeads.h>
#include <eigensym.h>
#include <assertext.h>

namespace EIGEN3EXT{
  
  /**
   * @class EigenSparseGenEigenSolver providing the interface for solving the
   * general eigen-value problem for the sparsematrix of the eigen3 library.
   * 
   * @see SparseGenEigenSolver
   * 
   * @note 
   * (1) only the lower part of K and M should be provided, otherwise there
   * will be error like "the coefffients in the matrix is not consistent" when
   * calling the arpack.
   * (2) the matrix K and M should be column compressed, which means that the 
   * function K.makeCompressed() and M.makeCompressed() shoulbe be called after 
   * initialization.
   */
  template <typename real>
  class EigenSparseGenEigenSolver{

  public:
	static bool solve (const Eigen::SparseMatrix<real> &K,
					   const Eigen::SparseMatrix<real> &M,
					   Eigen::Matrix<real,Eigen::Dynamic,Eigen::Dynamic> &eig_vec,
					   Eigen::Matrix<real,Eigen::Dynamic,1> &eig_val,
					   const int max_eig_num){

	  //check the parameters
	  assert(max_eig_num >0 && max_eig_num <= K.rows());
	  assert(K.nonZeros() > 0 && M.nonZeros() > 0);
	  assert_eq(K.rows(),M.rows());
	  assert_eq(K.cols(),M.cols());
	  assert_eq(K.rows(),K.cols());

	  //allocate memory
	  eig_vec.resize(K.rows(),max_eig_num);
	  eig_val.resize(max_eig_num);
  
	  //caculate 
	  const int k_nonzero_num = K.nonZeros();
	  const int k_rows = K.rows();
	  const double *k_data = K.valuePtr();
	  const int *k_rowind = K.innerIndexPtr();
	  const int *k_colptr = K.outerIndexPtr();

	  const int m_nonzero_num = M.nonZeros();
	  const double *m_data = M.valuePtr();
	  const int *m_rowind = M.innerIndexPtr();
	  const int *m_colptr = M.outerIndexPtr();

	  double *p_eig_val = &(eig_val[0]);
	  double *p_eig_vec = &(eig_vec(0,0));

	  bool succ = eigensym(k_data,k_nonzero_num,k_rowind,k_colptr,
	  					   m_data,m_nonzero_num,m_rowind,m_colptr,
	  					   k_rows,max_eig_num,
	  					   p_eig_val,p_eig_vec);

	  return  succ;
	}
	
  };

  typedef EigenSparseGenEigenSolver<double> EigenSparseGenEigenSolverD;
  typedef EigenSparseGenEigenSolver<float> EigenSparseGenEigenSolverF;
    
}//end of namespace

#endif /*_EIGENSPARSEGENEIGENSOLVER_H_*/
