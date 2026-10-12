#include <iostream>
#include "system.h"

using namespace std;

int main (int argc, char *argv[]){

  System SYS;

  SYS.initialize();
  SYS.initialize_properties();
  SYS.block_reset(0);

  // Loop over blocks
  for(int i=0; i<SYS.get_nbl(); i++){

    // Loop over simulation steps inside each block
    for(int j=0; j<SYS.get_nsteps(); j++){

      SYS.step();       // Perform one MC or MD step
      SYS.measure();    // Measure observables
    }

    SYS.averages(i+1);      // Block and progressive averages
    SYS.block_reset(i+1);   // Reset block accumulators
  }

  SYS.finalize();

  return 0;
}