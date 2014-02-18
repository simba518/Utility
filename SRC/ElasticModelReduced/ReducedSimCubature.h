#ifndef _REDUCEDSIMCUBATURE_H_
#define _REDUCEDSIMCUBATURE_H_

#include <stdlib.h>
#include <ElasticForceTetFullStVK.h>
#include "cubacode/GreedyCubop.h"
using namespace UTILITY;

namespace SIMULATOR{
  
  /**
   * @class ReducedSimCubature the cubature operator for the reduced simulation,
   * which compute the weights and sampled tetrahedrons.
   * 
   */
  class ReducedSimCubature:public GreedyCubop{
	
  public:
	ReducedSimCubature(const MatrixXd &B,pTetMesh_const tet_mesh):B(B),tet_mesh(tet_mesh){
	  assert(tet_mesh);
	  assert_eq(B.rows(), tet_mesh->nodes().size()*3);
	  tet_mesh->nodes(rest_shape);
	  fullStVKModel = pElasticForceTetFullStVK(new ElasticForceTetFullStVK(tet_mesh));
	}
	void run(const MatrixXd &trainingReduedDisp,
			 const MatrixXd &trainingFullForces,
			 double relErrTol,
			 int maxNumPoints,
			 int numCandsPerIter,
			 int itersPerFullNNLS,
			 int numSamplesPerSubtrain
			 ){

	  TrainingSet trainingSet(trainingReduedDisp.cols());
	  for (size_t i = 0; i < trainingSet.size(); ++i){
		trainingSet[i] = new VECTOR(trainingReduedDisp.rows());
		for (int j = 0; j < trainingSet[i]->size(); ++j)
		  (*trainingSet[i])(j) = trainingReduedDisp(j,i);
	  }
	  
	  VECTOR trainingForces(trainingFullForces.size());
	  memcpy(trainingForces.data(), &trainingFullForces(0,0), sizeof(double)*trainingFullForces.size());
	  
	  GreedyCubop::run(trainingSet, trainingForces, relErrTol, maxNumPoints, numCandsPerIter, itersPerFullNNLS, numSamplesPerSubtrain);
	  
	  for (size_t i = 0; i < trainingSet.size(); ++i){
		delete trainingSet[i];
	  }
	}
	const VectorXd &getWeights()const{
	  return weights;
	}
	const vector<int> &getSelectedTets()const{
	  return selectedTets;
	}
	
  protected:
	int numTotalPoints(){
	  return (NULL==tet_mesh) ? 0:tet_mesh->tets().size();
	}
	void evalPointForceDensity(int tet_id, VECTOR& q,VECTOR& gOut){
	  
	  assert_in(tet_id, 0, numTotalPoints()-1);
	  const VectorXd x = rest_shape+B*Map<VectorXd>(q.data(),q.size());

	  static mat3x4 f_tet;
	  fullStVKModel->force_tet(f_tet, tet_id, x);
	  VectorXd f(q.size());
	  f.setZero();
	  for (int j = 0;  j < 4; ++j){
		const int vi = tet_mesh->tets()[tet_id][j];
		f += B.block(3*vi,0,3,B.cols()).transpose()*((-1.0)*f_tet.col(j));
	  }

	  memcpy(gOut.data(), &f[0], sizeof(double)*f.size());
	}
	void handleCubature(vector<int>& selectedTets,VECTOR& weights,double relErr){
	  assert_eq(selectedTets.size(), weights.size());
	  assert_ge(relErr,0);
	  this->selectedTets = selectedTets;
	  this->weights = Map<VectorXd>(weights.data(), weights.size());
	  INFO_LOG("relative error of the cubature is: " <<relErr);
	}

  private:
	MatrixXd B;
	pTetMesh_const tet_mesh;
	VectorXd rest_shape;
	pElasticForceTetFullStVK fullStVKModel;
	VectorXd weights;
	vector<int> selectedTets;
  };
  
  typedef boost::shared_ptr<ReducedSimCubature> pReducedSimCubature;
  
}//end of namespace

#endif /*_REDUCEDSIMCUBATURE_H_*/
