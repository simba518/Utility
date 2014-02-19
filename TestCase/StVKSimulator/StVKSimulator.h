#ifndef _STVKSIMULATOR_H_
#define _STVKSIMULATOR_H_

#include <LagImpFullSim.h>
#include <FullStVKSimModel.h>
#include <Objmesh.h>
#include <Log.h>
using namespace UTILITY;

namespace SIMULATOR{
  
  /**
   * @class StVKSimulator a simulator using full stvk model.
   * 
   */
  class StVKSimulator{
	
  public:
	StVKSimulator(){

	  stvkModel = pFullStVKSimModel(new FullStVKSimModel());
	  simulator = pLagImpFullSim(new LagImpFullSim(stvkModel));
	  totalFrames = 200;
	  clearGravity();
	  saveRlstTo = "tempt/simulated_results";
	  steps = 1;
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
	void updateExtForces(VectorXd &fext)const;
	void buildFloor();
  
  private:
	pFullStVKSimModel stvkModel;
	pLagImpFullSim simulator;
	int totalFrames;
	int steps;
	vector<VectorXd> recorded_U;
	string saveRlstTo;
	double gravity[3];
	VectorXd f_ext;
	VectorXd u0;
	pObjmesh scene;
  };
  
  typedef boost::shared_ptr<StVKSimulator> pStVKSimulator;
  
}//end of namespace

#endif /*_STVKSIMULATOR_H_*/
