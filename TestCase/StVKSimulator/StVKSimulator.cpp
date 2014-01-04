#include "StVKSimulator.h"
#include <ConMatrixTools.h>
#include <MatrixIO.h>
#include <Timer.h>
using namespace UTILITY;
using namespace SIMULATOR;

// load data from init file: mesh, damping, time step, constraints, forces
void StVKSimulator::loadInitFile(const string filename){
  
  if(!simulator->init(filename)){
	ERROR_LOG("failed to initialze the simulator.");
  }
  
  UTILITY::JsonFilePaser jsonf;
  if (jsonf.open(filename)){
	string fixed_nodes_str;
	if(jsonf.readFilePath("fixed_nodes",fixed_nodes_str)){
	  vector<int> fixed_nodes;
	  UTILITY::loadVec(fixed_nodes_str, fixed_nodes, UTILITY::TEXT);
	  cout << "number of fixed nodes: " << fixed_nodes.size() << endl;
	  setFixedNodes(fixed_nodes);
	}

	if(!jsonf.read("save_results_to", saveRlstTo)){
	  saveRlstTo = "tempt/simulated_results";
	}
	if(!jsonf.read("T",totalFrames)){
	  totalFrames = 200;
	}

	vector<double> g;
	if(jsonf.read("gravity",g)){
	  assert_eq(g.size(),3);
	  gravity[0] = g[0];
	  gravity[1] = g[1];
	  gravity[2] = g[2];
	}else{
	  clearGravity();
	}

  }else{
	ERROR_LOG("failed to open the initfile: " << filename);
  }
}

void StVKSimulator::simulate(){

  bool succ = true;
  recorded_U.clear();
  recorded_U.reserve(totalFrames);
  cout << "total steps: " << totalFrames << endl;

  Timer timer;
  timer.start();
  for (int i = 0; i < totalFrames; ++i){

	cout << "step = " << i << endl;
	simulator->setExtForceForAllNodes(gravity[0], gravity[1], gravity[2]);
	const bool succ = simulator->forward();
	recorded_U.push_back(simulator->getU());

	ERROR_LOG_COND("simulation failed, step i = " << i, succ);
  }
  timer.stop("total simulation time is: ");
}

void StVKSimulator::save(){
  
  if(!EIGEN3EXT::write(saveRlstTo+".b",recorded_U)){
	ERROR_LOG("failed to save the results U to " << saveRlstTo<<".b");
  }else{
	cout << "success to save the results U to: " << saveRlstTo<<".b" << endl;
  }

  const pTetMesh_const tetmesh = stvkModel->getTetMesh();
  assert(tetmesh);
  if(!tetmesh->writeVTK(saveRlstTo, recorded_U)){
	ERROR_LOG("failed to save the results mesh to " << saveRlstTo<<"*.vtk");
  }else{
	cout << "success to save the results mesh to: " << saveRlstTo<<"*.vtk" << endl;
  }
  
}

void StVKSimulator::setFixedNodes(const vector<int> &fixednodes){

  VectorXd uc(fixednodes.size()*3);
  uc.setZero();
  simulator->setUc(uc);
  
  VecT trip_C;
  const int n = stvkModel->dimension()/3;
  UTILITY::computeConM(fixednodes, trip_C, n);
  simulator->setConM(trip_C);
}
