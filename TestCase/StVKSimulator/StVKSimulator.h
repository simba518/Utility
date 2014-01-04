#ifndef _STVKSIMULATOR_H_
#define _STVKSIMULATOR_H_

#include <LagImpFullSim.h>
#include <FullStVKSimModel.h>
#include <Log.h>

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
	}
	void loadInitFile(const string filename);
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
  
  private:
	pFullStVKSimModel stvkModel;
	pLagImpFullSim simulator;
	int totalFrames;
	vector<VectorXd> recorded_U;
	string saveRlstTo;
	double gravity[3];
  };
  
  typedef boost::shared_ptr<StVKSimulator> pStVKSimulator;
  
}//end of namespace

#endif /*_STVKSIMULATOR_H_*/
