#include "ReducedSolidSimulator.h"
using namespace SIMULATOR;

int main(int argc, char *argv[]){

  if (argc < 2){
	INFO_LOG("ussage: ReducedSimulator initfile");
	return 1;
  }

  ReducedSolidSimulator sim;
  sim.loadInitFile(argv[1]);
  sim.simulate();
  sim.save();
  return 0;
}
