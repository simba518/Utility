#ifndef _MASSMATRIX_H_
#define _MASSMATRIX_H_

#include <boost/shared_ptr.hpp>
#include <TetMesh.h>
#include <eigen3/Eigen/Sparse>

namespace UTILITY{

#define ADD_INFLUENCE(a,b)										\
  _M.coeffRef((int)tet[(int)a],(int)tet[(int)b])+=influence;	\
  _M.coeffRef((int)tet[(int)b],(int)tet[(int)a])+=influence;
  
  /// @class MassMatrix calculate the lumped mass matrix for tet mesh.
  /// @todo compute the umlumped mass matrix.
  class MassMatrix{
	
  public:
	void compute(Eigen::DiagonalMatrix<double,-1>&M,const TetMesh&mesh){

	  //assemble
	  _M.resize((int)mesh.nodes().size(),(int)mesh.nodes().size());
	  for(int i=0;i<(int)mesh.tets().size();i++)
	  	assembleMass(mesh,mesh.tets()[i],mesh.material()._rho[i]);
		
	  //accumulate entries
	  M.resize(mesh.nodes().size()*3);
	  M.setZero();
	  for(int k=0;k<_M.outerSize();++k)
	  	for(Eigen::SparseMatrix<double>::InnerIterator it(_M,k);it;++it)
	  	  M.diagonal().block<3,1>(it.row()*3,0)+=Vector3d::Constant(it.value());
	}
	static void lump(const Eigen::SparseMatrix<double> &fullM,
					 Eigen::DiagonalMatrix<double,-1>&M/*lumped M*/){
	  const int n3 = fullM.cols();
	  M.resize(n3);
	  M.setZero();
	  for (int i = 0; i < n3; ++i)
		M.diagonal()[i] = fullM.col(i).sum();
	}

  protected:
	void assembleMass(const TetMesh& mesh,const Vector4i& tet,const double& density){

	  double influence=tetrahedron(mesh.nodes()[tet[0]],mesh.nodes()[tet[1]],
								   mesh.nodes()[tet[2]],mesh.nodes()[tet[3]]).volume()/20.0f;
	  influence*=density;

	  ADD_INFLUENCE(0,0);
	  ADD_INFLUENCE(1,1);
	  ADD_INFLUENCE(2,2);
	  ADD_INFLUENCE(3,3);

	  ADD_INFLUENCE(1,0);
	  ADD_INFLUENCE(2,0);
	  ADD_INFLUENCE(3,0);

	  ADD_INFLUENCE(1,2);
	  ADD_INFLUENCE(2,3);
	  ADD_INFLUENCE(3,1);
	}

  private:
	Eigen::SparseMatrix<double> _M;
  };
  
  typedef boost::shared_ptr<MassMatrix> pMassMatrix;
  
}//end of namespace

#endif /* _MASSMATRIX_H_ */
