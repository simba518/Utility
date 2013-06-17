/**
 *\file eigensym.h
 *\brief examples for solving eigenvalue problems using arapack can be found at\n
 *http://www.inf-cv.uni-jena.de/proj/linal/eigen_8h-source.html\n
 *purpose:\n
 *solve the general eigen-problem: Ax = lambda*Bx with the n_eig largest eigen-values,
 *where A and B symetric and only the lower part of the matrix are stored.
 */
#ifndef _EIGENSYM_H_
#define _EIGENSYM_H_

/**
 *\fn bool eigensym( T *valA,int nnzA,  int *irowA,  int *pcolA,T *valB,int nnzB,  int *irowB,  int *pcolB,int n,int n_eig,T *vals,T *vecs)
 *
 *\brief  solve the general eigen-problem: Ax = lambda*Bx with the n_eig smallest eigen-values.
 *
 *\param[in] valA|valB pointer to an array that stores the nonzero elements of A and B.
 *\param[in] nnzA|nnzB Number of nonzero elements in A and B.
 *\param[in] irowA|irowB  pointer to an array that stores the row indices of the nonzeros in A and B.
 *\param[in] pcolA|pcolB pointer to an array of pointers to the beginning of each column of A (B) in valA (valB).
 *\param[in] n Dimension of the problem.
 *\param[in] n_eig number of eigenvectors to solve
 *
 *\param[out] vals  return the eigen values
 *\param[out] vecs  return the eigen vectors
 *\return success or not,and results is writed to vals and vecs
 *
 *\note
 *(1)the both the sparse matrix A and B should be symetric and only the lower part of the matrix are stored.
 *(2)the memory-space of vals and vecs should be allocated outside(before calling this function).\n
 *(3)onely float and double is supported, because of the limitation of arpack.\n
 *(4)number of the required eigen numbers should small than the dimension of the problem:\n
 *n_eig < n \n
 *\see eigensymTest.cpp
 */
namespace EIGEN3EXT{
  
  template <typename T>
  bool eigensym( const T *valA,int nnzA,  const int *irowA,  const int *pcolA,
				 const T *valB,int nnzB,  const int *irowB,  const int *pcolB,
				 int n,int n_eig,
				 T *vals,T *vecs);

}

#endif /* _EIGENSYM_H_ */
