#include "ReducedSimulator.h"
#include <JsonFilePaser.h>
#include <Log.h>
using namespace UTILITY;
using namespace SIMULATOR;

bool ReducedSimulator::init(const string filename){

  JsonFilePaser jsonf;
  bool succ = false;
  if (jsonf.open(filename)){
	succ = true;
	this->reset();
	jsonf.read("alpha_k",alpha_k,0.0);
	jsonf.read("alpha_m",alpha_m,0.0);
	jsonf.read("h",h,1.0);
  }else{
	ERROR_LOG("failed to open the initfile: " << filename);
  }
  return succ;
}

void ReducedImpLogConSimulator::setConGroups(const vector<set<int> > &groups){
  
  
}

void ReducedImpLogConSimulator::setUc(const VectorXd &uc){

  
}

void ReducedImpLogConSimulator::removeAllCon(){
  

}

bool ReducedImpLogConSimulator::forward(){
  

}
