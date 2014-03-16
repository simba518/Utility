#ifndef _COMPUTEMASSMAT_H_
#define _COMPUTEMASSMAT_H_

#include <TetMesh.h>
#include <eigen3/Eigen/Sparse>
#include <CASADITools.h>
using namespace UTILITY;
using CasADi::SX;
using CasADi::SXMatrix;

namespace ELASTIC_OPT{

#define ADD_INFLUENCE_AD(a,b,M)								\
  M.elem((int)tet[(int)a],(int)tet[(int)b])+=influence;		\
  M.elem((int)tet[(int)b],(int)tet[(int)a])+=influence;
  
  /**
   * @class compute lumped mass matrix using casadi.
   * 
   */
  class ComputeMassMat{

  public:
	void compute(SXMatrix &M,const TetMesh&mesh){
	  
	  vector<SX> density(mesh.material()._rho.size());
	  for (int i = 0;i<mesh.material()._rho.size(); ++i){
		density[i] = mesh.material()._rho[i];
	  }
	  compute(M,mesh, density);
	}
	void compute(SXMatrix &M,const TetMesh&mesh, const vector<SX> &density){
	  
	  assert_eq(density.size(),mesh.tets().size());
	  this->_density = density;

	  _M.setZero();
	  computeCompactM(_M,mesh);

	  const int n = mesh.nodes().size();
	  M.resize(n*3,n*3);
	  M.setZero();

	  for (int i = 0; i < _M.size1(); ++i){
		for (int j = 0; j < _M.size2(); ++j){
		  if( _M.hasNZ(i,j) ){
			M.elem(i*3+0,i*3+0) += _M.elem(i,j);
			M.elem(i*3+1,i*3+1) += _M.elem(i,j);
			M.elem(i*3+2,i*3+2) += _M.elem(i,j);
		  }
		}
	  }
	}

  protected:
	void assembleMass(const TetMesh& mesh,const Vector4i& tet,const SX& density,SXMatrix &M)const{

	  SX influence=tetrahedron(mesh.nodes()[tet[0]],mesh.nodes()[tet[1]],
								   mesh.nodes()[tet[2]],mesh.nodes()[tet[3]]).volume()/20.0f;
	  influence*=density;

	  ADD_INFLUENCE_AD(0,0,M);
	  ADD_INFLUENCE_AD(1,1,M);
	  ADD_INFLUENCE_AD(2,2,M);
	  ADD_INFLUENCE_AD(3,3,M);

	  ADD_INFLUENCE_AD(1,0,M);
	  ADD_INFLUENCE_AD(2,0,M);
	  ADD_INFLUENCE_AD(3,0,M);

	  ADD_INFLUENCE_AD(1,2,M);
	  ADD_INFLUENCE_AD(2,3,M);
	  ADD_INFLUENCE_AD(3,1,M);
	}
	void computeCompactM(SXMatrix &M,const TetMesh& mesh)const{

	  const int n = (int)mesh.nodes().size();
	  M.resize(n,n);
	  for(int i=0;i<(int)mesh.tets().size();i++){
	  	assembleMass(mesh,mesh.tets()[i],_density[i],M);
	  }
	}

  private:
	SXMatrix _M; // compact mass matrix.
	vector<SX> _density;
  };

  typedef boost::shared_ptr<ComputeMassMat> pComputeMassMat;
  
}//end of namespace

#endif /*_COMPUTEMASSMAT_H_*/
