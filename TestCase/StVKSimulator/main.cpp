#include "StVKSimulator.h"

int main(int argc, char *argv[]){

  if (argc < 2){
	INFO_LOG("ussage: StVKSimulator initfile");
	INFO_LOG("example: ./StVKSimulator Data/dino/simulate.ini");
	return 1;
  }

  SIMULATOR::StVKSimulator sim;
  sim.loadInitFile(argv[1]);
  sim.simulate();
  sim.save();
  return 0;
}
