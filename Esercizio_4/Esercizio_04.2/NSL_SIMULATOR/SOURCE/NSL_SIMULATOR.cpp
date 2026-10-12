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

  // Initial low-entropy configuration: t* = 0
  SYS.write_XYZ(0);

  for(int i=0; i<SYS.get_nbl(); i++){

    for(int j=0; j<SYS.get_nsteps(); j++){

      SYS.step();
      SYS.measure();

      global_step++;

      // Save selected configurations during the evolution
      if(global_step == 100    ||   // t* = 0.1
         global_step == 200    ||   // t* = 0.2
         global_step == 300    ||   // t* = 0.3
         global_step == 500    ||   // t* = 0.5
         global_step == 1000   ||   // t* = 1
         global_step == 5000   ||   // t* = 5
         global_step == 20000  ||   // t* = 20
         global_step == 50000  ||   // t* = 50
         global_step == 100000){    // t* = 100

        SYS.write_XYZ(global_step);
      }
    }

    SYS.averages(i+1);
    SYS.block_reset(i+1);
  }

  SYS.finalize();

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
