#ifndef _REDUCEDSOLIDSIMULATOR_H_
#define _REDUCEDSOLIDSIMULATOR_H_

#include <ReducedSimulator.h>
#include <Objmesh.h>
#include <Log.h>
using namespace UTILITY;

namespace SIMULATOR{
  
  /**
   * @class ReducedSolidSimulator a example for simulation of solid using reduced method.
   * 
   */
  class ReducedSolidSimulator{
	
  public:
	ReducedSolidSimulator(){

	  totalFrames = 200;
	  clearGravity();
	  saveRlstTo = "tempt/simulated_results";
	  steps = 1;
	  tetmesh = pTetMesh(new TetMesh());
	}
	void loadInitFile(const string filename);
	void initCollision();
	void setGravity(const double g[3]){
	  gravity[0] = g[0];
	  gravity[1] = g[1];
	  gravity[2] = g[2];
	}
	void clearGravity(){
	  gravity[0] = 0.0f;
	  gravity[1] = 0.0f;
	  gravity[2] = 0.0f;	  
	}
	void simulate();
	void save();

  protected:
	void setFixedNodes(const vector<int> &fixednodes);
	void updateExtForces(VectorXd &fext);
	void buildFloor();
  
  private:
	pTetMesh tetmesh;
	pCubaturedElasticModel stvkModel;
	pReducedSimulator simulator;
	int totalFrames;
	int steps;
	vector<VectorXd> recorded_U;
	string saveRlstTo;
	double gravity[3];
	VectorXd f_ext;
	VectorXd u0, u, x;
	pObjmesh scene;
  };
  
  typedef boost::shared_ptr<ReducedSolidSimulator> pReducedSolidSimulator;
  
}//end of namespace

#endif /*_REDUCEDSOLIDSIMULATOR_H_*/
