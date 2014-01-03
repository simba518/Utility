#include <LagImpFullSim.h>
#include <FullStVKSimModel.h>
#include <Log.h>
using namespace SIMULATOR;

class FullSimulator{
  
public:
  FullSimulator(){

	tetmesh = pTetMesh(new TetMesh());
	stvkModel = pFullStVKSimModel(new FullStVKSimModel(tetmesh));
	simulator = pLagImpFullSim(new LagImpFullSim(stvkModel));
	totalFrames = 200;
  }
  void loadInitFile(const string filename){
	// load data from init file: mesh, damping, time step, constraints, forces
  }
  void simulate(){
	
  }
  void save(){
	
  }

protected:
  void setFixedNodes(const vector<int> &fixednodes){
	
  }
  void setGravity(const double g[3]){
	
  }
  
private:
  pTetMesh tetmesh;
  pFullStVKSimModel stvkModel;
  pLagImpFullSim simulator;

  int totalFrames;
  vector<VectorXd> U;
  string saveRlstTo;
};

int main(int argc, char *argv[]){

  FullSimulator sim;
  assert_ge(argc,2);
  sim.loadInitFile(argv[1]);
  sim.simulate();
  sim.save();
  return 0;
}
