#include "StVKSimulator.h"
using namespace SIMULATOR;

int main(int argc, char *argv[]){

  if (argc < 2){
	INFO_LOG("ussage: StVKSimulator initfile");
	return 1;
  }

  StVKSimulator sim;
  sim.loadInitFile(argv[1]);
  sim.simulate();
  sim.save();
  return 0;
}
