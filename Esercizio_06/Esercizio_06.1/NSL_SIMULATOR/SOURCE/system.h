/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __System__
#define __System__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <armadillo>
#include <stdlib.h>
#include "particle.h"
#include "random.h"

using namespace std;
using namespace arma;

class System {

private:
  const int _ndim = 3;
  bool _restart;
  int _sim_type;
  int _npart;
  int _nblocks;
  int _nsteps;
  int _neqblocks;
  int _neqsteps;
  int _nattempts;
  int _naccepted;
  double _temp, _beta;
  double _rho, _volume;
  double _r_cut;
  double _delta;
  double _J, _H;
  string _start_state;
  string _run_mode;
  string _run_tag;

  vec _side;
  vec _halfside;
  Random _rnd;
  field<Particle> _particle;
  vec _fx, _fy, _fz;

  // Proprieta' misurate
  int _nprop;
  bool _measure_penergy, _measure_kenergy, _measure_tenergy;
  bool _measure_temp, _measure_pressure, _measure_gofr;
  bool _measure_magnet, _measure_cv, _measure_chi;
  bool _measure_pofv;

  int _index_penergy, _index_kenergy, _index_tenergy;
  int _index_temp, _index_pressure, _index_gofr;
  int _index_magnet, _index_cv, _index_chi;
  int _index_pofv;

  int _n_bins;
  int _n_bins_v;
  double _bin_size;
  double _bin_size_v;
  double _vtail, _ptail;

  vec _block_av;
  vec _global_av;
  vec _global_av2;
  vec _average;
  vec _measurement;

  string number_tag(double value, int precision) const;
  string build_run_tag() const;
  string output_file(const string &stem) const;

public:
  // Getter utili nel main
  int get_nbl();
  int get_nsteps();
  int get_neqbl();
  int get_neqsteps();
  int get_sim_type();
  double get_temp();
  double get_field();
  string get_start_state();
  string get_run_mode();
  string get_run_tag();

  // Riconfigurazione usata da RUN_MODE ALL_EQ
  void configure_ising_run(
    int sim_type,
    double temp,
    double field,
    const string &start_state
);

  // Inizializzazione e I/O
  void initialize();
  void initialize_properties();
  void finalize();
  void write_configuration();
  void write_XYZ(int nconf);
  void read_configuration();
  void initialize_velocities();

  // Dinamica
  void step();
  void move(int part);
  bool metro(int part);
  void Verlet();

  // Misure
  void block_reset(int blk);
  void measure();
  void averages(int blk);
  double error(double acc, double acc2, int blk);
  double ising_energy_per_spin();
  double ising_magnetization_per_spin();

  // Funzioni di supporto
  double pbc(double position, int i);
  int pbc(int i);
  double Force(int i, int dim);
  double Boltzmann(int i, bool xnew);
};

#endif // __System__
