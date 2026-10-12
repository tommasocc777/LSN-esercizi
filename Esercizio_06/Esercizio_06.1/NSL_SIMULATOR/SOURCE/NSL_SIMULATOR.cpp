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
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

#include "system.h"

using namespace std;

// Converte la temperatura nel formato usato dal notebook:
// 0.5 -> T0p5, 1.0 -> T1p0, 2.0 -> T2p0.
string temperature_tag(double T){
  ostringstream ss;
  ss << fixed << setprecision(1) << T;
  string tag = ss.str();
  replace(tag.begin(), tag.end(), '.', 'p');
  return "T" + tag;
}

// Salva il transiente di equilibrazione.
// L'energia viene registrata ad ogni MC step e, in parallelo,
// viene calcolata una media locale su NEQSTEPS step.
// Salva il transiente di equilibrazione.
//
// Energia e magnetizzazione per spin vengono registrate
// ad ogni MC step.
//
// In parallelo vengono calcolate le corrispondenti
// medie locali su NEQSTEPS MC steps.
//
// I dati istantanei verranno utilizzati nel notebook
// anche per costruire C/N e chi/N dalle fluttuazioni.
void run_equilibration(System &SYS,
                       const string &energy_file,
                       const string &blocks_file){

  ofstream eq_inst(energy_file);
  ofstream eq_block(blocks_file);

  if(!eq_inst.is_open() || !eq_block.is_open()){

    cerr
      << "PROBLEMA: impossibile aprire i file di equilibrazione"
      << endl;

    exit(EXIT_FAILURE);
  }

  // Intestazione dei file
  eq_inst
    << "# MC_STEP:    U/N:    M/N:"
    << endl;

  eq_block
    << "# BLOCK:    MC_STEPS:    BLOCK_U/N:    BLOCK_M/N:"
    << endl;

  int mc_step = 0;


  // ============================================================
  // STATO INIZIALE
  // ============================================================

  // Vengono salvate energia e magnetizzazione
  // prima della prima mossa Monte Carlo.
  eq_inst
    << setw(12) << mc_step
    << setw(16) << SYS.ising_energy_per_spin()
    << setw(16) << SYS.ising_magnetization_per_spin()
    << endl;


  // ============================================================
  // EQUILIBRAZIONE
  // ============================================================

  for(int iblk=0; iblk<SYS.get_neqbl(); iblk++){

    double block_energy = 0.0;
    double block_magnetization = 0.0;


    for(int istep=0; istep<SYS.get_neqsteps(); istep++){

      // Viene eseguito un Monte Carlo step
      SYS.step();

      mc_step++;


      // Vengono misurate energia e magnetizzazione per spin
      const double energy =
        SYS.ising_energy_per_spin();

      const double magnetization =
        SYS.ising_magnetization_per_spin();


      // Vengono salvati i valori istantanei
      eq_inst
        << setw(12) << mc_step
        << setw(16) << energy
        << setw(16) << magnetization
        << endl;


      // Vengono aggiornati gli accumulatori locali
      block_energy += energy;
      block_magnetization += magnetization;
    }


    // Vengono calcolate le medie locali
    block_energy /=
      double(SYS.get_neqsteps());

    block_magnetization /=
      double(SYS.get_neqsteps());


    // Vengono salvate le medie locali
    eq_block
      << setw(12) << iblk + 1
      << setw(14) << mc_step
      << setw(16) << block_energy
      << setw(16) << block_magnetization
      << endl;
  }


  eq_inst.close();
  eq_block.close();
}

// ============================================================================
// EQUILIBRAZIONE SENZA OUTPUT
// ============================================================================

// Vengono eseguiti i MC steps necessari all'equilibrazione.
// I dati di questa fase non vengono salvati, poiche' lo studio
// dell'equilibrazione e' gia' stato effettuato nel notebook.
void equilibrate(System &SYS){

  for(int iblk=0; iblk<SYS.get_neqbl(); iblk++){

    for(int istep=0; istep<SYS.get_neqsteps(); istep++){

      SYS.step();
    }
  }
}

// ============================================================================
// PRODUZIONE
// ============================================================================

// Dopo l'equilibrazione vengono inizializzati gli accumulatori
// e viene effettuato il data blocking delle osservabili.
void run_production(System &SYS){

  // Gli accumulatori vengono inizializzati solamente
  // dopo il completamento dell'equilibrazione.
  SYS.initialize_properties();

  SYS.block_reset(0);


  for(int iblk=0; iblk<SYS.get_nbl(); iblk++){

    for(int istep=0; istep<SYS.get_nsteps(); istep++){

      SYS.step();

      SYS.measure();
    }


    SYS.averages(
        iblk + 1
    );

    SYS.block_reset(
        iblk + 1
    );
  }
}



int main(){

  System SYS;

  // Viene inizializzato il sistema leggendo i file di input.
  SYS.initialize();

  const string mode = SYS.get_run_mode();


  // ========================================================================
  // RUN_MODE ALL_EQ
  // ========================================================================
  //
  // Vengono eseguite automaticamente le simulazioni utilizzate
  // per lo studio preliminare dell'equilibrazione:
  //
  //   T = 0.5, 1.0, 2.0
  //
  // per:
  //
  //   - Metropolis
  //   - Gibbs
  //   - configurazione ordinata
  //   - configurazione casuale
  //
  // Il campo magnetico viene posto a zero.
  //
  if(mode == "all_eq"){

    const vector<int> simulation_types = {
      2,  // Metropolis
      3   // Gibbs
    };

    const vector<int> temperatures_10 = {
      5,   // T = 0.5
      10,  // T = 1.0
      20   // T = 2.0
    };

    const vector<string> states = {
      "ordered",
      "random"
    };


    for(int sim_type : simulation_types){

      const string method =
          (sim_type == 2)
          ? "metropolis"
          : "gibbs";


      for(int temp10 : temperatures_10){

        const double T =
            double(temp10) / 10.0;

        const string T_tag =
            "T"
            + to_string(temp10 / 10)
            + "p"
            + to_string(temp10 % 10);


        for(const string &state : states){

          cout
            << "\n----------------------------------------"
            << endl;

          cout
            << method
            << " | T = "
            << T
            << " | "
            << state
            << endl;


          // Viene configurata la nuova simulazione.
          SYS.configure_ising_run(
              sim_type,
              T,
              0.0,
              state
          );


          const string energy_file =
              "../OUTPUT/"
              "equilibration_energy_"
              + method
              + "_"
              + T_tag
              + "_"
              + state
              + ".dat";


          const string blocks_file =
              "../OUTPUT/"
              "equilibration_blocks_"
              + method
              + "_"
              + T_tag
              + "_"
              + state
              + ".dat";


          run_equilibration(
              SYS,
              energy_file,
              blocks_file
          );
        }
      }
    }


    SYS.finalize();

    cout
      << "\nTutte le simulazioni di equilibrazione "
      << "sono state completate."
      << endl;

    return 0;
  }


  // ========================================================================
  // RUN_MODE ALL_PROD
  // ========================================================================
  //
  // Vengono eseguite automaticamente tutte le simulazioni necessarie
  // per costruire le curve termodinamiche:
  //
  //   T = 0.5, 0.6, ..., 2.0
  //
  // Per ogni temperatura:
  //
  //   - Metropolis
  //   - Gibbs
  //
  // e:
  //
  //   h = 0.00  -> U/N, C/N, chi/N
  //   h = 0.02  -> M/N
  //
  // Ogni simulazione parte da una configurazione casuale.
  //
  if(mode == "all_prod"){

    const vector<int> simulation_types = {
      2,  // Metropolis
      3   // Gibbs
    };

    const vector<double> fields = {
      0.00,
      0.02
    };


    // Vengono considerate 16 temperature:
    //
    // 0.5, 0.6, ..., 2.0
    //
    for(int itemp=0; itemp<=15; itemp++){

      const double T =
          0.5
          + 0.1 * double(itemp);


      for(int sim_type : simulation_types){

        const string method =
            (sim_type == 2)
            ? "Metropolis"
            : "Gibbs";


        for(double h : fields){

          cout
            << "\n----------------------------------------"
            << endl;

          cout
            << method
            << " | T = "
            << fixed
            << setprecision(1)
            << T
            << " | h = "
            << setprecision(2)
            << h
            << endl;


          // Viene inizializzata una nuova configurazione casuale.
          SYS.configure_ising_run(
              sim_type,
              T,
              h,
              "random"
          );


          // ================================================================
          // EQUILIBRAZIONE
          // ================================================================

          cout
            << "  Equilibrazione..."
            << endl;

          equilibrate(
              SYS
          );


          // ================================================================
          // PRODUZIONE E DATA BLOCKING
          // ================================================================

          cout
            << "  Produzione..."
            << endl;

          run_production(
              SYS
          );


          cout
            << "  Completata."
            << endl;
        }
      }
    }


    SYS.finalize();

    cout
      << "\nTutte le simulazioni di produzione "
      << "sono state completate."
      << endl;

    return 0;
  }


  // ========================================================================
  // RUN_MODE EQ_ONLY
  // ========================================================================
  //
  // Viene effettuata solamente la fase di equilibrazione.
  // Energia e magnetizzazione vengono salvate durante il transiente.
  //
  if(mode == "eq_only"){

    run_equilibration(
        SYS,
        "../OUTPUT/equilibration_energy.dat",
        "../OUTPUT/equilibration_blocks.dat"
    );

    SYS.finalize();

    return 0;
  }


  // ========================================================================
  // RUN_MODE SINGLE
  // ========================================================================
  //
  // Viene effettuata una singola simulazione completa:
  //
  //   1. equilibrazione;
  //   2. produzione;
  //   3. data blocking.
  //
  if(mode == "single"){

    run_equilibration(
        SYS,
        "../OUTPUT/equilibration_energy.dat",
        "../OUTPUT/equilibration_blocks.dat"
    );


    run_production(
        SYS
    );


    SYS.finalize();

    return 0;
  }


  // ========================================================================
  // CONTROLLO DEL RUN_MODE
  // ========================================================================

  cerr
    << "PROBLEMA: RUN_MODE non riconosciuto."
    << endl;

  cerr
    << "Utilizzare SINGLE, EQ_ONLY, ALL_EQ oppure ALL_PROD."
    << endl;


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
