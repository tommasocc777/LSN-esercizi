/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cctype>
#include "system.h"

using namespace std;
using namespace arma;

// ============================================================================
// FUNZIONI DI SUPPORTO PER I NOMI DEI FILE
// ============================================================================

string System :: number_tag(double value, int precision) const {
  ostringstream ss;
  ss << fixed << setprecision(precision) << value;
  string out = ss.str();

  replace(out.begin(), out.end(), '.', 'p');
  replace(out.begin(), out.end(), '-', 'm');

  return out;
}

string System :: build_run_tag() const {
  string method;

  if(_sim_type == 2) method = "metropolis";
  else if(_sim_type == 3) method = "gibbs";
  else method = "simulation";

  string state = _restart ? "restart" : _start_state;

  return method
       + "_T" + number_tag(_temp, 1)
       + "_" + state
       + "_h" + number_tag(_H, 3);
}

string System :: output_file(const string &stem) const {

  // Nella produzione automatica ogni simulazione deve avere
  // un file distinto, identificato da algoritmo, temperatura,
  // configurazione iniziale e campo magnetico.
  if(_run_mode == "all_prod"){
    return "../OUTPUT/"
         + stem
         + "_"
         + _run_tag
         + ".dat";
  }

  // Per SINGLE, EQ_ONLY e ALL_EQ viene mantenuto
  // il comportamento utilizzato finora.
  return "../OUTPUT/" + stem + ".dat";
}

// ============================================================================
// DINAMICA
// ============================================================================

void System :: step(){
  if(_sim_type == 0) {
    this->Verlet();
  } else {
    // Un MC step corrisponde a N tentativi di aggiornamento.
    for(int i=0; i<_npart; i++) {
      this->move(int(_rnd.Rannyu()*_npart));
    }
  }

  _nattempts += _npart;
}

void System :: move(int i){

  if(_sim_type == 3){
    // -----------------------------------------------------------------------
    // GIBBS / HEAT-BATH PER ISING 1D
    // -----------------------------------------------------------------------
    // Il nuovo spin viene estratto direttamente dalla distribuzione
    // condizionata fissati i due vicini.

    const int left  = this->pbc(i-1);
    const int right = this->pbc(i+1);

    const double h_loc =
        _J * (_particle(left).getspin() + _particle(right).getspin()) + _H;

    const double p_up = 1.0 / (1.0 + exp(-2.0 * _beta * h_loc));

    if(_rnd.Rannyu() < p_up) _particle(i).setspin(+1);
    else                     _particle(i).setspin(-1);

    // Nel Gibbs non esiste una fase di rifiuto della proposta.
    _naccepted++;
  }
  else {
    // -----------------------------------------------------------------------
    // METROPOLIS
    // -----------------------------------------------------------------------
    if(_sim_type == 1){
      vec shift(_ndim);

      for(int j=0; j<_ndim; j++) {
        shift(j) = _rnd.Rannyu(-1.0,1.0) * _delta;
      }

      _particle(i).translate(shift, _side);

      if(this->metro(i)){
        _particle(i).acceptmove();
        _naccepted++;
      } else {
        _particle(i).moveback();
      }
    }
    else {
      // Ising 1D con single-spin flip.
      if(this->metro(i)){
        _particle(i).flip();
        _naccepted++;
      }
    }
  }
}

bool System :: metro(int i){
  double delta_E;

  if(_sim_type == 1) {
    delta_E = this->Boltzmann(i,true) - this->Boltzmann(i,false);
  }
  else {
    delta_E = 2.0 * _particle(i).getspin()
            * (_J * (_particle(this->pbc(i-1)).getspin()
                   + _particle(this->pbc(i+1)).getspin())
               + _H);
  }

  const double acceptance = min(1.0, exp(-_beta * delta_E));

  return (_rnd.Rannyu() < acceptance);
}

void System :: Verlet(){
  double xnew, ynew, znew;

  for(int i=0; i<_npart; i++){
    _fx(i) = this->Force(i,0);
    _fy(i) = this->Force(i,1);
    _fz(i) = this->Force(i,2);
  }

  for(int i=0; i<_npart; i++){
    xnew = this->pbc(2.0 * _particle(i).getposition(0,true)
                   - _particle(i).getposition(0,false)
                   + _fx(i) * pow(_delta,2), 0);

    ynew = this->pbc(2.0 * _particle(i).getposition(1,true)
                   - _particle(i).getposition(1,false)
                   + _fy(i) * pow(_delta,2), 1);

    znew = this->pbc(2.0 * _particle(i).getposition(2,true)
                   - _particle(i).getposition(2,false)
                   + _fz(i) * pow(_delta,2), 2);

    _particle(i).setvelocity(0,
      this->pbc(xnew - _particle(i).getposition(0,false), 0)/(2.0*_delta));
    _particle(i).setvelocity(1,
      this->pbc(ynew - _particle(i).getposition(1,false), 1)/(2.0*_delta));
    _particle(i).setvelocity(2,
      this->pbc(znew - _particle(i).getposition(2,false), 2)/(2.0*_delta));

    _particle(i).acceptmove();
    _particle(i).setposition(0, xnew);
    _particle(i).setposition(1, ynew);
    _particle(i).setposition(2, znew);
  }

  _naccepted += _npart;
}

// ============================================================================
// ENERGIA E FUNZIONI DI SUPPORTO
// ============================================================================

double System :: ising_energy_per_spin(){
  double energy = 0.0;

  for(int i=0; i<_npart; i++){
    const double s_i = double(_particle(i).getspin());
    const double s_j = double(_particle(this->pbc(i+1)).getspin());

    energy += -_J * s_i * s_j - 0.5 * _H * (s_i + s_j);
  }

  return energy / double(_npart);
}

// Restituisce la magnetizzazione totale divisa per il numero di spin.
// La funzione viene utilizzata anche durante la fase di equilibrazione.
double System :: ising_magnetization_per_spin(){

  double magnetization = 0.0;

  for(int i=0; i<_npart; i++){
    magnetization += double(_particle(i).getspin());
  }

  return magnetization / double(_npart);
}
double System :: Boltzmann(int i, bool xnew){
  double energy_i = 0.0;
  double dx, dy, dz, dr;

  for(int j=0; j<_npart; j++){
    if(j != i){
      dx = this->pbc(_particle(i).getposition(0,xnew)
                   - _particle(j).getposition(0,true), 0);
      dy = this->pbc(_particle(i).getposition(1,xnew)
                   - _particle(j).getposition(1,true), 1);
      dz = this->pbc(_particle(i).getposition(2,xnew)
                   - _particle(j).getposition(2,true), 2);

      dr = sqrt(dx*dx + dy*dy + dz*dz);

      if(dr < _r_cut) {
        energy_i += 1.0/pow(dr,12) - 1.0/pow(dr,6);
      }
    }
  }

  return 4.0 * energy_i;
}

double System :: Force(int i, int dim){
  double f = 0.0;
  double dr;
  vec distance(_ndim);

  for(int j=0; j<_npart; j++){
    if(i != j){
      distance(0) = this->pbc(_particle(i).getposition(0,true)
                            - _particle(j).getposition(0,true), 0);
      distance(1) = this->pbc(_particle(i).getposition(1,true)
                            - _particle(j).getposition(1,true), 1);
      distance(2) = this->pbc(_particle(i).getposition(2,true)
                            - _particle(j).getposition(2,true), 2);

      dr = sqrt(dot(distance,distance));

      if(dr < _r_cut){
        f += distance(dim) * (48.0/pow(dr,14) - 24.0/pow(dr,8));
      }
    }
  }

  return f;
}

double System :: pbc(double position, int i){
  return position - _side(i) * rint(position/_side(i));
}

int System :: pbc(int i){
  if(i >= _npart) i -= _npart;
  else if(i < 0) i += _npart;
  return i;
}

// ============================================================================
// INIZIALIZZAZIONE
// ============================================================================

void System :: initialize(){

  // Valori di default.
  _restart = false;
  _sim_type = 2;
  _npart = 0;
  _nblocks = 0;
  _nsteps = 0;
  _neqblocks = 0;
  _neqsteps = 0;
  _nattempts = 0;
  _naccepted = 0;
  _temp = 1.0;
  _beta = 1.0;
  _rho = 1.0;
  _volume = 1.0;
  _r_cut = 0.0;
  _delta = 0.0;
  _J = 1.0;
  _H = 0.0;
  _start_state = "random";
  _run_mode = "single";

  // Inizializzazione del generatore di numeri casuali.
  int p1, p2;
  ifstream Primes("../INPUT/Primes");

  if(!Primes.is_open()){
    cerr << "PROBLEMA: impossibile aprire ../INPUT/Primes" << endl;
    exit(EXIT_FAILURE);
  }

  Primes >> p1 >> p2;
  Primes.close();

  int seed[4];
  ifstream Seed("../INPUT/seed.in");

  if(!Seed.is_open()){
    cerr << "PROBLEMA: impossibile aprire ../INPUT/seed.in" << endl;
    exit(EXIT_FAILURE);
  }

  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  Seed.close();

  _rnd.SetRandom(seed,p1,p2);

  // Lettura dei parametri della simulazione.
  ifstream input("../INPUT/input.dat");

  if(!input.is_open()){
    cerr << "PROBLEMA: impossibile aprire ../INPUT/input.dat" << endl;
    exit(EXIT_FAILURE);
  }

  string property;
  double delta;

  while(input >> property){

    if(property == "RUN_MODE"){
      input >> _run_mode;

      transform(_run_mode.begin(), _run_mode.end(),
                _run_mode.begin(),
                [](unsigned char c){ return char(tolower(c)); });
    }
    else if(property == "SIMULATION_TYPE"){
      input >> _sim_type;

      if(_sim_type > 1){
        input >> _J >> _H;
      }

      if(_sim_type < 0 || _sim_type > 3){
        cerr << "PROBLEMA: simulation type non valido" << endl;
        exit(EXIT_FAILURE);
      }
    }
    else if(property == "RESTART"){
      input >> _restart;
    }
    else if(property == "START_STATE"){
      input >> _start_state;

      transform(_start_state.begin(), _start_state.end(),
                _start_state.begin(),
                [](unsigned char c){ return char(tolower(c)); });
    }
    else if(property == "TEMP"){
      input >> _temp;
      _beta = 1.0/_temp;
    }
    else if(property == "NPART"){
      input >> _npart;

      _fx.resize(_npart);
      _fy.resize(_npart);
      _fz.resize(_npart);
      _particle.set_size(_npart);

      // Particle::initialize() imposta inizialmente tutti gli spin a +1.
      for(int i=0; i<_npart; i++) {
        _particle(i).initialize();
      }
    }
    else if(property == "RHO"){
      input >> _rho;
      _volume = _npart/_rho;

      _side.resize(_ndim);
      _halfside.resize(_ndim);

      const double side = pow(_volume,1.0/3.0);

      for(int i=0; i<_ndim; i++) _side(i) = side;
      _halfside = 0.5*_side;
    }
    else if(property == "R_CUT"){
      input >> _r_cut;
    }
    else if(property == "DELTA"){
      input >> delta;
      _delta = delta;
    }
    else if(property == "NEQBLOCKS"){
      input >> _neqblocks;
    }
    else if(property == "NEQSTEPS"){
      input >> _neqsteps;
    }
    else if(property == "NBLOCKS"){
      input >> _nblocks;
    }
    else if(property == "NSTEPS"){
      input >> _nsteps;
    }
    else if(property == "ENDINPUT"){
      break;
    }
    else {
      cerr << "PROBLEMA: input sconosciuto -> " << property << endl;
      exit(EXIT_FAILURE);
    }
  }

  input.close();

  if(_npart <= 0){
    cerr << "PROBLEMA: NPART deve essere positivo" << endl;
    exit(EXIT_FAILURE);
  }

  // Configurazione iniziale del modello di Ising.
  if(_sim_type > 1){

    if(_restart){
      this->read_configuration();
    }
    else if(_start_state == "ordered"){
      for(int i=0; i<_npart; i++) _particle(i).setspin(+1);
    }
    else if(_start_state == "random"){
      for(int i=0; i<_npart; i++){
        if(_rnd.Rannyu() < 0.5) _particle(i).setspin(+1);
        else                    _particle(i).setspin(-1);
      }
    }
    else {
      cerr << "PROBLEMA: START_STATE deve essere ordered oppure random"
           << endl;
      exit(EXIT_FAILURE);
    }
  }
  else {
    this->read_configuration();
  }

  if(_sim_type == 0) this->initialize_velocities();

  _run_tag = this->build_run_tag();

  // File di acceptance specifico per questa run.
  ofstream couta(this->output_file("acceptance"));
  couta << "#   N_BLOCK:  ACCEPTANCE:" << endl;
  couta.close();

  // Log specifico per questa run.
  ofstream coutf(this->output_file("output"));

  if(_sim_type == 2) coutf << "ISING 1D MONTE CARLO (METROPOLIS)" << endl;
  else if(_sim_type == 3) coutf << "ISING 1D MONTE CARLO (GIBBS)" << endl;
  else if(_sim_type == 1) coutf << "LJ MONTE CARLO (NVT)" << endl;
  else coutf << "LJ MOLECULAR DYNAMICS (NVE)" << endl;

  coutf << "RUN_MODE= " << _run_mode << endl;
  coutf << "RUN_TAG= " << _run_tag << endl;
  coutf << "TEMPERATURE= " << _temp << endl;
  coutf << "J= " << _J << endl;
  coutf << "H= " << _H << endl;
  coutf << "START_STATE= " << (_restart ? "restart" : _start_state) << endl;
  coutf << "NPART= " << _npart << endl;
  coutf << "NEQBLOCKS= " << _neqblocks << endl;
  coutf << "NEQSTEPS= " << _neqsteps << endl;
  coutf << "NBLOCKS= " << _nblocks << endl;
  coutf << "NSTEPS= " << _nsteps << endl;
  coutf << "System initialized!" << endl;
  coutf.close();
}

void System :: read_configuration(){

  // Per il modello di Ising il file config.spin viene letto soltanto
  // quando RESTART=1.
  if(_sim_type > 1){

    if(!_restart) return;

    ifstream cinf("../INPUT/CONFIG/config.spin");

    if(!cinf.is_open()){
      cerr << "PROBLEMA: impossibile aprire ../INPUT/CONFIG/config.spin"
           << endl;
      exit(EXIT_FAILURE);
    }

    for(int i=0; i<_npart; i++){
      int spin;

      if(!(cinf >> spin)){
        cerr << "PROBLEMA: config.spin contiene meno di "
             << _npart << " spin" << endl;
        exit(EXIT_FAILURE);
      }

      if(spin != +1 && spin != -1){
        cerr << "PROBLEMA: spin non valido in config.spin -> "
             << spin << endl;
        exit(EXIT_FAILURE);
      }

      _particle(i).setspin(spin);
    }

    cinf.close();
    return;
  }

  // Lettura della configurazione per Lennard-Jones.
  ifstream cinf("../INPUT/CONFIG/config.xyz");

  if(!cinf.is_open()){
    cerr << "PROBLEMA: impossibile aprire ../INPUT/CONFIG/config.xyz"
         << endl;
    exit(EXIT_FAILURE);
  }

  string comment, particle;
  double x, y, z;
  int ncoord;

  cinf >> ncoord;

  if(ncoord != _npart){
    cerr << "PROBLEMA: numero di coordinate non compatibile con NPART"
         << endl;
    exit(EXIT_FAILURE);
  }

  cinf >> comment;

  for(int i=0; i<_npart; i++){
    cinf >> particle >> x >> y >> z;

    _particle(i).setposition(0, this->pbc(_side(0)*x,0));
    _particle(i).setposition(1, this->pbc(_side(1)*y,1));
    _particle(i).setposition(2, this->pbc(_side(2)*z,2));
    _particle(i).acceptmove();
  }

  cinf.close();
}

void System :: initialize_velocities(){
  double xold, yold, zold;

  if(_restart){
    ifstream cinf("../INPUT/CONFIG/conf-1.xyz");

    if(!cinf.is_open()){
      cerr << "PROBLEMA: impossibile aprire ../INPUT/CONFIG/conf-1.xyz"
           << endl;
      exit(EXIT_FAILURE);
    }

    string comment, particle;
    int ncoord;

    cinf >> ncoord;

    if(ncoord != _npart){
      cerr << "PROBLEMA: numero di coordinate non compatibile con NPART"
           << endl;
      exit(EXIT_FAILURE);
    }

    cinf >> comment;

    for(int i=0; i<_npart; i++){
      cinf >> particle >> xold >> yold >> zold;
      _particle(i).setpositold(0,this->pbc(_side(0)*xold,0));
      _particle(i).setpositold(1,this->pbc(_side(1)*yold,1));
      _particle(i).setpositold(2,this->pbc(_side(2)*zold,2));
    }

    cinf.close();
  }
  else {
    vec vx(_npart), vy(_npart), vz(_npart);
    vec sumv(_ndim, fill::zeros);

    for(int i=0; i<_npart; i++){
      vx(i) = _rnd.Gauss(0.0,sqrt(_temp));
      vy(i) = _rnd.Gauss(0.0,sqrt(_temp));
      vz(i) = _rnd.Gauss(0.0,sqrt(_temp));

      sumv(0) += vx(i);
      sumv(1) += vy(i);
      sumv(2) += vz(i);
    }

    sumv /= double(_npart);

    double sumv2 = 0.0;

    for(int i=0; i<_npart; i++){
      vx(i) -= sumv(0);
      vy(i) -= sumv(1);
      vz(i) -= sumv(2);

      sumv2 += vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i);
    }

    sumv2 /= double(_npart);
    const double scalef = sqrt(3.0*_temp/sumv2);

    for(int i=0; i<_npart; i++){
      _particle(i).setvelocity(0,vx(i)*scalef);
      _particle(i).setvelocity(1,vy(i)*scalef);
      _particle(i).setvelocity(2,vz(i)*scalef);
    }

    for(int i=0; i<_npart; i++){
      xold = this->pbc(_particle(i).getposition(0,true)
                       - _particle(i).getvelocity(0)*_delta,0);
      yold = this->pbc(_particle(i).getposition(1,true)
                       - _particle(i).getvelocity(1)*_delta,1);
      zold = this->pbc(_particle(i).getposition(2,true)
                       - _particle(i).getvelocity(2)*_delta,2);

      _particle(i).setpositold(0,xold);
      _particle(i).setpositold(1,yold);
      _particle(i).setpositold(2,zold);
    }
  }
}

// ============================================================================
// PROPRIETA' E DATA BLOCKING
// ============================================================================

void System :: initialize_properties(){

  string property;
  int index_property = 0;
  _nprop = 0;

  _measure_penergy = false;
  _measure_kenergy = false;
  _measure_tenergy = false;
  _measure_temp = false;
  _measure_pressure = false;
  _measure_gofr = false;
  _measure_magnet = false;
  _measure_cv = false;
  _measure_chi = false;
  _measure_pofv = false;

  ifstream input("../INPUT/properties.dat");

  if(!input.is_open()){
    cerr << "PROBLEMA: impossibile aprire ../INPUT/properties.dat" << endl;
    exit(EXIT_FAILURE);
  }

  while(input >> property){

    if(property == "POTENTIAL_ENERGY"){
      ofstream out(this->output_file("potential_energy"));
      out << "#     BLOCK:  ACTUAL_PE:     PE_AVE:      ERROR:" << endl;
      out.close();

      _nprop++;
      _index_penergy = index_property++;
      _measure_penergy = true;
      _vtail = 0.0;
    }
    else if(property == "KINETIC_ENERGY"){
      ofstream out(this->output_file("kinetic_energy"));
      out << "#     BLOCK:   ACTUAL_KE:    KE_AVE:      ERROR:" << endl;
      out.close();

      _nprop++;
      _index_kenergy = index_property++;
      _measure_kenergy = true;
    }
    else if(property == "TOTAL_ENERGY"){
      ofstream out(this->output_file("total_energy"));
      out << "#     BLOCK:   ACTUAL_TE:    TE_AVE:      ERROR:" << endl;
      out.close();

      _nprop++;
      _index_tenergy = index_property++;
      _measure_tenergy = true;
    }
    else if(property == "TEMPERATURE"){
      ofstream out(this->output_file("temperature"));
      out << "#     BLOCK:   ACTUAL_T:     T_AVE:       ERROR:" << endl;
      out.close();

      _nprop++;
      _index_temp = index_property++;
      _measure_temp = true;
    }
    else if(property == "PRESSURE"){
      ofstream out(this->output_file("pressure"));
      out << "#     BLOCK:   ACTUAL_P:     P_AVE:       ERROR:" << endl;
      out.close();

      _nprop++;
      _index_pressure = index_property++;
      _measure_pressure = true;
      _ptail = 0.0;
    }
    else if(property == "GOFR"){
      input >> _n_bins;
      _nprop += _n_bins;
      _bin_size = _halfside.min()/double(_n_bins);
      _measure_gofr = true;
      _index_gofr = index_property;
      index_property += _n_bins;
    }
    else if(property == "POFV"){
      if(_sim_type > 0){
        cerr << "PROBLEMA: POFV non definita per simulazioni MC" << endl;
        exit(EXIT_FAILURE);
      }

      input >> _n_bins_v;
      _nprop += _n_bins_v;
      _bin_size_v = 4.0*_temp/double(_n_bins_v);
      _measure_pofv = true;
      _index_pofv = index_property;
      index_property += _n_bins_v;
    }
    else if(property == "MAGNETIZATION"){
      ofstream out(this->output_file("magnetization"));
      out << "#     BLOCK:   ACTUAL_M:     M_AVE:       ERROR:" << endl;
      out.close();

      _nprop++;
      _index_magnet = index_property++;
      _measure_magnet = true;
    }
    else if(property == "SPECIFIC_HEAT"){
      ofstream out(this->output_file("specific_heat"));
      out << "#     BLOCK:   ACTUAL_CV:    CV_AVE:      ERROR:" << endl;
      out.close();

      _nprop++;
      _index_cv = index_property++;
      _measure_cv = true;
    }
    else if(property == "SUSCEPTIBILITY"){
      ofstream out(this->output_file("susceptibility"));
      out << "#     BLOCK:   ACTUAL_X:     X_AVE:       ERROR:" << endl;
      out.close();

      _nprop++;
      _index_chi = index_property++;
      _measure_chi = true;
    }
    else if(property == "ENDPROPERTIES"){
      break;
    }
    else {
      cerr << "PROBLEMA: proprieta' sconosciuta -> " << property << endl;
      exit(EXIT_FAILURE);
    }
  }

  input.close();

  if(_measure_cv && !_measure_tenergy){
    cerr << "PROBLEMA: SPECIFIC_HEAT richiede anche TOTAL_ENERGY "
         << "in properties.dat" << endl;
    exit(EXIT_FAILURE);
  }

  _measurement.resize(_nprop);
  _average.resize(_nprop);
  _block_av.resize(_nprop);
  _global_av.resize(_nprop);
  _global_av2.resize(_nprop);

  _measurement.zeros();
  _average.zeros();
  _block_av.zeros();
  _global_av.zeros();
  _global_av2.zeros();

  // I contatori vengono azzerati dopo l'equilibrazione.
  _nattempts = 0;
  _naccepted = 0;

  ofstream log(this->output_file("output"), ios::app);
  log << "Reading properties completed!" << endl;
  log.close();
}

void System :: measure(){
  _measurement.zeros();

  double penergy_temp = 0.0;
  double kenergy_temp = 0.0;
  double virial = 0.0;

  // -------------------------------------------------------------------------
  // Lennard-Jones
  // -------------------------------------------------------------------------
  if(_measure_penergy || _measure_pressure || _measure_gofr){
    vec distance(_ndim);

    for(int i=0; i<_npart-1; i++){
      for(int j=i+1; j<_npart; j++){
        distance(0) = this->pbc(_particle(i).getposition(0,true)
                              - _particle(j).getposition(0,true),0);
        distance(1) = this->pbc(_particle(i).getposition(1,true)
                              - _particle(j).getposition(1,true),1);
        distance(2) = this->pbc(_particle(i).getposition(2,true)
                              - _particle(j).getposition(2,true),2);

        const double dr = sqrt(dot(distance,distance));

        if(dr < _r_cut){
          if(_measure_penergy)
            penergy_temp += 1.0/pow(dr,12) - 1.0/pow(dr,6);

          if(_measure_pressure)
            virial += 1.0/pow(dr,12) - 0.5/pow(dr,6);
        }
      }
    }
  }

  if(_measure_penergy){
    penergy_temp = _vtail + 4.0*penergy_temp/double(_npart);
    _measurement(_index_penergy) = penergy_temp;
  }

  if(_measure_kenergy){
    for(int i=0; i<_npart; i++){
      kenergy_temp += 0.5 * dot(_particle(i).getvelocity(),
                                _particle(i).getvelocity());
    }

    kenergy_temp /= double(_npart);
    _measurement(_index_kenergy) = kenergy_temp;
  }

  // -------------------------------------------------------------------------
  // Energia
  // -------------------------------------------------------------------------
  double u = 0.0;

  if(_sim_type > 1 && (_measure_tenergy || _measure_cv)){
    u = this->ising_energy_per_spin();
  }

  if(_measure_tenergy){
    if(_sim_type < 2)
      _measurement(_index_tenergy) = kenergy_temp + penergy_temp;
    else
      _measurement(_index_tenergy) = u;
  }

  if(_measure_temp && _measure_kenergy)
    _measurement(_index_temp) = (2.0/3.0)*kenergy_temp;

  if(_measure_pressure)
    _measurement(_index_pressure) =
      _rho*(2.0/3.0)*kenergy_temp
      + (_ptail*_npart + 48.0*virial/3.0)/_volume;

  // -------------------------------------------------------------------------
  // Magnetizzazione e suscettibilita' per Ising
  // -------------------------------------------------------------------------
  if(_sim_type > 1 && (_measure_magnet || _measure_chi)){
    double m_tot = 0.0;

    for(int i=0; i<_npart; i++)
      m_tot += double(_particle(i).getspin());

    if(_measure_magnet)
      _measurement(_index_magnet) = m_tot/double(_npart);

    // A h=0 la suscettibilita' per spin richiesta dall'esercizio e'
    // chi/N = beta * <M^2> / N.
    if(_measure_chi)
      _measurement(_index_chi) = _beta*m_tot*m_tot/double(_npart);
  }

  // Per C/N viene accumulato u^2; la combinazione non lineare viene
  // effettuata alla fine del blocco.
  if(_sim_type > 1 && _measure_cv)
    _measurement(_index_cv) = u*u;

  _block_av += _measurement;
}

void System :: averages(int blk){

  ofstream coutf;
  double average, sum_average, sum_ave2;

  _average = _block_av/double(_nsteps);

  // La capacita' termica viene costruita dalle fluttuazioni all'interno
  // del singolo blocco:
  // C/N = beta^2 * N * (<u^2> - <u>^2).
  if(_measure_cv){
    const double u_mean  = _average(_index_tenergy);
    const double u2_mean = _average(_index_cv);

    _average(_index_cv) =
      _beta*_beta*double(_npart)*(u2_mean - u_mean*u_mean);
  }

  _global_av  += _average;
  _global_av2 += _average % _average;

  auto write_property = [&](bool enabled, int index, const string &stem){
    if(!enabled) return;

    coutf.open(this->output_file(stem), ios::app);

    average = _average(index);
    sum_average = _global_av(index);
    sum_ave2 = _global_av2(index);

    coutf << setw(12) << blk
          << setw(16) << average
          << setw(16) << sum_average/double(blk)
          << setw(16) << this->error(sum_average,sum_ave2,blk)
          << endl;

    coutf.close();
  };

  write_property(_measure_penergy,  _index_penergy, "potential_energy");
  write_property(_measure_kenergy,  _index_kenergy, "kinetic_energy");
  write_property(_measure_tenergy,  _index_tenergy, "total_energy");
  write_property(_measure_temp,     _index_temp,    "temperature");
  write_property(_measure_pressure, _index_pressure,"pressure");
  write_property(_measure_magnet,   _index_magnet,  "magnetization");
  write_property(_measure_cv,       _index_cv,      "specific_heat");
  write_property(_measure_chi,      _index_chi,     "susceptibility");

  // Rapporto di accettazione progressivo.
  double fraction = 0.0;

  if(_nattempts > 0)
    fraction = double(_naccepted)/double(_nattempts);

  coutf.open(this->output_file("acceptance"),ios::app);
  coutf << setw(12) << blk
        << setw(16) << fraction
        << endl;
  coutf.close();
}

void System :: block_reset(int blk){
  if(blk > 0){
    ofstream log(this->output_file("output"), ios::app);
    log << "Block completed: " << blk << endl;
    log.close();
  }

  _block_av.zeros();
}

double System :: error(double acc, double acc2, int blk){

  if(blk <= 1){
    return 0.0;
  }

  return sqrt(
      fabs(
          acc2 / double(blk)
          - pow(
              acc / double(blk),
              2
          )
      )
      / double(blk - 1)
  );
}
// ============================================================================
// OUTPUT FINALE
// ============================================================================

void System :: write_configuration(){

  if(_sim_type > 1){
    ofstream out("../OUTPUT/CONFIG/config.spin");

    for(int i=0; i<_npart; i++)
      out << _particle(i).getspin() << " ";

    out << endl;
    out.close();
    return;
  }

  ofstream out("../OUTPUT/CONFIG/config.xyz");

  out << _npart << endl;
  out << "#Comment!" << endl;

  for(int i=0; i<_npart; i++){
    out << "LJ  "
        << setprecision(17) << _particle(i).getposition(0,true)/_side(0) << "   "
        << setprecision(17) << _particle(i).getposition(1,true)/_side(1) << "   "
        << setprecision(17) << _particle(i).getposition(2,true)/_side(2) << endl;
  }

  out.close();
}

void System :: write_XYZ(int nconf){
  ofstream out("../OUTPUT/CONFIG/config_" + to_string(nconf) + ".xyz");

  if(!out.is_open()) return;

  out << _npart << endl;
  out << "#Comment!" << endl;

  for(int i=0; i<_npart; i++){
    out << "LJ  "
        << setw(16) << _particle(i).getposition(0,true)
        << setw(16) << _particle(i).getposition(1,true)
        << setw(16) << _particle(i).getposition(2,true)
        << endl;
  }

  out.close();
}

void System :: finalize(){
  this->write_configuration();
  _rnd.SaveSeed();

  ofstream log(this->output_file("output"),ios::app);
  log << "Simulation completed!" << endl;
  log.close();
}

// ============================================================================
// RICONFIGURAZIONE PER LE 12 RUN DI EQUILIBRAZIONE
// ============================================================================

void System :: configure_ising_run(
    int sim_type,
    double temp,
    double field,
    const string &start_state
){

  if(sim_type != 2 && sim_type != 3){

    cerr
      << "PROBLEMA: configure_ising_run richiede "
      << "SIMULATION_TYPE 2 o 3"
      << endl;

    exit(EXIT_FAILURE);
  }


  // Vengono aggiornati i parametri della simulazione
  _sim_type = sim_type;

  _temp = temp;
  _beta = 1.0 / _temp;

  _H = field;

  _restart = false;
  _start_state = start_state;


  transform(
      _start_state.begin(),
      _start_state.end(),
      _start_state.begin(),
      [](unsigned char c){
        return char(tolower(c));
      }
  );


  // ============================================================
  // CONFIGURAZIONE INIZIALE
  // ============================================================

  if(_start_state == "ordered"){

    for(int i=0; i<_npart; i++){
      _particle(i).setspin(+1);
    }

  }
  else if(_start_state == "random"){

    for(int i=0; i<_npart; i++){

      if(_rnd.Rannyu() < 0.5)
        _particle(i).setspin(+1);

      else
        _particle(i).setspin(-1);
    }

  }
  else {

    cerr
      << "PROBLEMA: START_STATE deve essere "
      << "ORDERED oppure RANDOM"
      << endl;

    exit(EXIT_FAILURE);
  }


  // Vengono azzerati i contatori
  _nattempts = 0;
  _naccepted = 0;


  // Viene costruito il nome identificativo della nuova run
  _run_tag = this->build_run_tag();


  // ============================================================
  // FILE SPECIFICI DELLA RUN AUTOMATICA
  // ============================================================

  if(_run_mode == "all_prod"){

    // File del rapporto di accettazione
    ofstream acceptance(
        this->output_file("acceptance")
    );

    acceptance
      << "#   N_BLOCK:  ACCEPTANCE:"
      << endl;

    acceptance.close();


    // File di log
    ofstream log(
        this->output_file("output")
    );

    log
      << "RUN_TAG= "
      << _run_tag
      << endl;

    log
      << "TEMPERATURE= "
      << _temp
      << endl;

    log
      << "H= "
      << _H
      << endl;

    log.close();
  }
}

// ============================================================================
// GETTER
// ============================================================================

int System :: get_nbl(){ return _nblocks; }
int System :: get_nsteps(){ return _nsteps; }
int System :: get_neqbl(){ return _neqblocks; }
int System :: get_neqsteps(){ return _neqsteps; }
int System :: get_sim_type(){ return _sim_type; }
double System :: get_temp(){ return _temp; }
double System :: get_field(){ return _H; }
string System :: get_start_state(){ return _restart ? "restart" : _start_state; }
string System :: get_run_mode(){ return _run_mode; }
string System :: get_run_tag(){ return _run_tag; }

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
