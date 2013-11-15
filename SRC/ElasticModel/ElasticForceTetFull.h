#ifndef _ELASTICFORCETETFULL_H_
#define _ELASTICFORCETETFULL_H_

#include <boost/shared_ptr.hpp>
#include <eigen3/Eigen/Sparse>
#include <DefGradTet.h>

namespace UTILITY{

  struct TetDF{ // dF/dX for one tetrahedron.
    Matrix3d df[4][4];
  };
  typedef Eigen::Matrix<double,3,4> mat3x4;
  typedef std::vector<mat3x4,Eigen::aligned_allocator<mat3x4> > VecMat3x4;
  
  /**
   * @class ElasticForceTetFull base class for full elastic forces of
   * tetrahedron mesh.
   * 
   */
  class ElasticForceTetFull{
	
  public:
	ElasticForceTetFull(pTetMesh_const vol_mesh):_vol_mesh(vol_mesh){
	  this->prepare();
	}
	virtual void force(const VectorXd &X, VectorXd &out);
	virtual const SparseMatrix<double> &K(const VectorXd &X);
	virtual void Kdx(const VectorXd &dx, const VectorXd &X, VectorXd &out);

  protected:
	virtual void prepare();
	virtual void computeTetForces(VecMat3x4 &tet_forces, const VectorXd &X)=0;
	virtual void computeTetForceDerivX(std::vector<TetDF>&df,const VectorXd &X)=0;
	virtual void computeTetForceDerivXdX(VecMat3x4 &tet_kdx,const VectorXd &dx,const VectorXd &X)=0;
	void initKEntry(SparseMatrix<double> &K,std::vector<int> &entries)const;
	
  protected:
	pTetMesh_const _vol_mesh;
	DefGradTet _def_grad;
	std::vector<double> _volume;
	VecMat3x4 _tet_forces;
	std::vector<TetDF> _tet_k;
	VecMat3x4 _tet_kdx;
	SparseMatrix<double> _Kx; // K(x)
	std::vector<int> _entries; 
  };
  
  typedef boost::shared_ptr<ElasticForceTetFull> pElasticForceTetFull;
  
}//end of namespace

#endif /*_ELASTICFORCETETFULL_H_*/
