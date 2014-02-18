#include "ReducedElasticModel.h"
#include <MassMatrix.h>
#include <JsonFilePaser.h>
using namespace UTILITY;
using namespace SIMULATOR;

bool ReducedElasticModel::init(const string filename){
  
  JsonFilePaser jsonf;
  bool succ = false;
  if (jsonf.open(filename)){
  	succ = jsonf.readMatFile("subspace_basis",B);
  }else{
  	ERROR_LOG("failed to open the initfile: " << filename);
  }
  return succ;
}

bool CubaturedElasticModel::init(const string filename){

  bool succ = ReducedElasticModel::init(filename);
  if (!succ){
  	return false;
  }

  JsonFilePaser jsonf;
  if (jsonf.open(filename)){
  	succ = jsonf.readVecFile("cubature_weights",weights);
  	succ &= jsonf.readVecFile("cubature_points",sampledTets,UTILITY::TEXT);

  	UTILITY::MassMatrix mass;
  	Eigen::DiagonalMatrix<double,-1> diag_M;
  	assert(ElasticForceTetFull::_vol_mesh);
  	mass.compute(diag_M, *ElasticForceTetFull::_vol_mesh);
  	assert_eq(diag_M.rows(),B.rows());
  	M = B.transpose()*diag_M*B;
  }else{
  	ERROR_LOG("failed to open the initfile: " << filename);
  }
  return succ;
}

// f = sum wi*Bi^t*fi(q), for i in S.
bool CubaturedElasticModel::evaluateF(const VectorXd &reduced_u,VectorXd &f){
  
  assert_eq(weights.size(), sampledTets.size());
  assert_eq(rest_shape.size(), B.rows());
  assert_eq(reduced_u.size(), B.cols());
  const VectorXd x = rest_shape+B*reduced_u;

  f.resize(reducedDim());
  f.setZero();

  static mat3x4 f_tet;
  for (size_t i = 0; i < sampledTets.size(); ++i){
	const int tet_id = sampledTets[i];
    force_tet(f_tet, tet_id, x);
	for (int j = 0;  j < 4; ++j){
	  const int vi = _vol_mesh->tets()[tet_id][j];
	  f += B.block(3*vi,0,3,B.cols()).transpose()*((-1.0f*weights[i])*f_tet.col(j));
	}
  }
  return true;
}

// K = sum wi*Bi^t*Ki(q)*Bi, for i in S.
bool CubaturedElasticModel::evaluateK(const VectorXd &reduced_u,MatrixXd &K){
  
  assert_eq(weights.size(), sampledTets.size());
  assert_eq(rest_shape.size(), B.rows());
  assert_eq(reduced_u.size(), B.cols()); 
  const VectorXd x = rest_shape+B*reduced_u;

  K.resize(reducedDim(),reducedDim());
  K.setZero();

  TetDF K_tet;
  for (size_t i = 0; i < sampledTets.size(); ++i){
    forceDerivX_tet(K_tet, sampledTets[i], x);
	for (int j = 0; j < 4; ++j){
	  const int vj = _vol_mesh->tets()[i][j];
	  const MatrixXd &Bj = B.block(3*vj,0,3,B.cols());
	  for (int k = 0; k < 4; ++k){
		const int vk = _vol_mesh->tets()[i][k];
		const MatrixXd &Bk = B.block(3*vk,0,3,B.cols());
		const Matrix3d Ajk = (-1.0f*weights[i])*K_tet.df[j][k];
		K += Bj.transpose()*Ajk*Bk;
	  }
	}
  }
  return true;
}
