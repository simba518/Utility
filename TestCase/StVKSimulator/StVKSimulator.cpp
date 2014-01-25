#include "StVKSimulator.h"
#include <ConMatrixTools.h>
#include <MatrixIO.h>
#include <Timer.h>
#include <CCD_DPF.h>
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

	if(!jsonf.readFilePath("save_results_to", saveRlstTo,false)){
	  saveRlstTo = "tempt/simulated_results";
	}
	if(!jsonf.read("T",totalFrames)){
	  totalFrames = 200;
	}
	if(!jsonf.read("steps",steps)){
	  steps = 1;
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

	jsonf.read("gravity_start", gravity_start);
	jsonf.read("gravity_stop", gravity_stop);
	jsonf.readVecFile("u0",u0);

	double response_scalor = 1.0f;
	if (jsonf.read("response_scalor",response_scalor)){
	  CCD_DPF::getInstance()->setReponseScalor(response_scalor);
	}

  }else{
	ERROR_LOG("failed to open the initfile: " << filename);
  }
  buildFloor();
  initCollision();
}

void StVKSimulator::initCollision(){

  CCD_DPF::getInstance()->removeAllCollisionObjects();
  if (scene && scene->getVertsNum() > 0){
	CCD_DPF::getInstance()->addObjectV(scene.get(),scene->getVerts(),scene->getFaces());
  }
  pTetMesh_const tetmesh = stvkModel->getTetMesh();
  if (tetmesh && tetmesh->nodes().size()>0){
	CCD_DPF::getInstance()->addObject(tetmesh.get(),tetmesh->nodes(),tetmesh->surface());
  }
  CCD_DPF::getInstance()->prepare();
}

void StVKSimulator::simulate(){

  bool succ = true;
  recorded_U.clear();
  recorded_U.reserve(totalFrames);
  cout << "total steps: " << totalFrames << endl;
  cout << "gravity start: " << gravity_start << endl;
  cout << "gravity stop: " << gravity_stop << endl;

  Timer timer;
  timer.start();
  if (u0.size() > 0){
	assert_eq(u0.size(), stvkModel->dimension());
	simulator->setU0(u0);
  }

  for (int i = 0; i < totalFrames; ++i){

	cout << "step = " << i << endl;

	bool succ = false;
	for (int k = 0; k < steps; ++k){
	  updateExtForces(f_ext);
	  simulator->setExtForce(f_ext);
	  succ = simulator->forward();
	}

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

  if (fixednodes.size() <= 0){
	simulator->removeAllCon();
	return;
  }

  VectorXd uc(fixednodes.size()*3);
  uc.setZero();
  simulator->setUc(uc);
  
  VecT trip_C;
  const int n = stvkModel->dimension()/3;
  UTILITY::computeConM(fixednodes, trip_C, n);
  simulator->setConM(trip_C);
}

void StVKSimulator::updateExtForces(VectorXd &f_ext)const{
  
  pTetMesh_const tetmesh = stvkModel->getTetMesh();
  assert(tetmesh);
  f_ext.resize(tetmesh->nodes().size()*3);
  for (int k = 0; k < f_ext.size()/3; ++k){
	f_ext[k*3+0] = gravity[0];
	f_ext[k*3+1] = gravity[1];
	f_ext[k*3+2] = gravity[2];
  }
  
  VectorXd x;
  tetmesh->nodes(x);
  x += simulator->getU();
  CCD_DPF::getInstance()->updateV(tetmesh.get(),x);
  CCD_DPF::getInstance()->checkCollision();
  CCD_DPF::getInstance()->computeForce();
  CCD_DPF::getInstance()->addCollisionForce(tetmesh.get(),f_ext);
}

void StVKSimulator::buildFloor(){
  
  const double y = -2.0;
  const int w_2 = 3.0f;
  VectorXd nodes(12);
  nodes << w_2,y,w_2,  -w_2,y,w_2,  -w_2,y,-w_2,  w_2,y,-w_2;
  VectorXi faces(6);
  faces << 0,2,1,  0,3,2;
  scene = pObjmesh(new Objmesh(nodes, faces));
  scene->write("./tempt/floor.obj");
}
