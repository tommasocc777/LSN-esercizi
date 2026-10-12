/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include "system.h"

using namespace std;

int main (int argc, char *argv[]){

  System SYS;

  SYS.initialize();
  SYS.initialize_properties();
  SYS.block_reset(0);

  int global_step = 0;

  int nblocks = SYS.get_nbl();
  int nsteps = SYS.get_nsteps();

  int forward_steps = nblocks * nsteps;

  // Save the initial low-entropy configuration
  SYS.write_XYZ(0);


  // ============================================================
  // FORWARD EVOLUTION
  // ============================================================

  cout << "Forward evolution..." << endl;

  for(int i=0; i<nblocks; i++){

    for(int j=0; j<nsteps; j++){

      SYS.step();
      SYS.measure();

      global_step++;
    }

    SYS.averages(i+1);
    SYS.block_reset(i+1);
  }


  // Save the state reached at the end of the forward trajectory
  SYS.write_XYZ(forward_steps);


  // ============================================================
  // TIME REVERSAL
  // ============================================================

  cout << "Reversing time at step "
       << forward_steps << "..." << endl;

  SYS.reverse_time();


  // ============================================================
  // BACKWARD EVOLUTION
  // ============================================================

  cout << "Backward evolution..." << endl;

  for(int i=0; i<nblocks; i++){

    for(int j=0; j<nsteps; j++){

      SYS.step();
      SYS.measure();

      global_step++;
    }

    // Continue numbering blocks after the forward evolution
    int block = nblocks + i + 1;

    SYS.averages(block);
    SYS.block_reset(block);
  }


  // Configuration obtained after returning backward
  SYS.write_XYZ(2 * forward_steps);

  SYS.finalize();

  cout << "Simulation completed." << endl;

  return 0;
}
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
