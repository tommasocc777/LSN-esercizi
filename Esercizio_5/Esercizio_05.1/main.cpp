#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <limits>

#include "random.h"

using namespace std;


// ============================================================
// Struttura per rappresentare una posizione nello spazio 3D
// ============================================================

struct Position {

    double x;
    double y;
    double z;

};


// ============================================================
// Struttura contenente i risultati finali di una simulazione
// ============================================================

struct SimulationResult {

    double mean;
    double error;
    double acceptance;

};


// ============================================================
// Calcolo della distanza dall'origine
// ============================================================

double Radius(const Position& p) {

    return sqrt(
        p.x * p.x +
        p.y * p.y +
        p.z * p.z
    );
}


// ============================================================
// Logaritmo della densità di probabilità
//
// Nell'algoritmo di Metropolis non è necessario conoscere
// le costanti di normalizzazione, perché si cancellano nel
// rapporto tra probabilità nuova e probabilità vecchia.
// ============================================================

double LogProbability(const Position& p,
                      const string& orbital) {

    const double r = Radius(p);


    // --------------------------------------------------------
    // Stato fondamentale 1s
    //
    // |Psi_100|^2 proporzionale a exp(-2r)
    //
    // log(p) = -2r + costante
    // --------------------------------------------------------

    if (orbital == "1s") {

        return -2.0 * r;
    }


    // --------------------------------------------------------
    // Stato eccitato 2p
    //
    // |Psi_210|^2 proporzionale a z^2 * exp(-r)
    //
    // log(p) = 2 log|z| - r + costante
    // --------------------------------------------------------

    if (orbital == "2p") {

        // Sul piano nodale z = 0 la densità di probabilità
        // è nulla
        if (fabs(p.z) < 1.0e-14) {

            return -numeric_limits<double>::infinity();
        }

        return 2.0 * log(fabs(p.z)) - r;
    }


    cerr << "Orbitale non riconosciuto." << endl;

    exit(1);
}


// ============================================================
// Generazione di una nuova configurazione proposta
//
// proposal = "uniform"
//     spostamento uniforme in [-step, step]
//
// proposal = "gauss"
//     spostamento gaussiano N(0, step^2)
//
// Le tre coordinate vengono generate indipendentemente.
// ============================================================

Position ProposeMove(const Position& position,
                     double step,
                     const string& proposal,
                     Random& rnd) {

    Position trial;


    // --------------------------------------------------------
    // Transition probability uniforme
    // --------------------------------------------------------

    if (proposal == "uniform") {

        trial.x =
            position.x + rnd.Rannyu(-step, step);

        trial.y =
            position.y + rnd.Rannyu(-step, step);

        trial.z =
            position.z + rnd.Rannyu(-step, step);


        return trial;
    }


    // --------------------------------------------------------
    // Transition probability gaussiana
    //
    // Ogni coordinata viene spostata secondo una normale
    // con media zero e deviazione standard step.
    // --------------------------------------------------------

    if (proposal == "gauss") {

        trial.x =
            position.x + rnd.Gauss(0.0, step);

        trial.y =
            position.y + rnd.Gauss(0.0, step);

        trial.z =
            position.z + rnd.Gauss(0.0, step);


        return trial;
    }


    cerr << "Transition probability non riconosciuta."
         << endl;

    exit(1);
}


// ============================================================
// Singolo passo dell'algoritmo di Metropolis
//
// Sia la proposal uniforme sia quella gaussiana utilizzata
// sono simmetriche:
//
// T(r'|r) = T(r|r')
//
// quindi la probabilità di accettazione è:
//
// A = min[1, p(r') / p(r)]
//
// La funzione restituisce true se la mossa viene accettata
// e false se viene rifiutata.
// ============================================================

bool MetropolisStep(Position& position,
                    double step,
                    const string& orbital,
                    const string& proposal,
                    Random& rnd) {

    // Generazione della nuova configurazione proposta
    const Position trial =
        ProposeMove(
            position,
            step,
            proposal,
            rnd
        );


    // --------------------------------------------------------
    // Logaritmo delle probabilità
    // --------------------------------------------------------

    const double log_p_old =
        LogProbability(position, orbital);

    const double log_p_new =
        LogProbability(trial, orbital);


    // --------------------------------------------------------
    // Logaritmo del rapporto tra le probabilità
    // --------------------------------------------------------

    const double log_ratio =
        log_p_new - log_p_old;


    // --------------------------------------------------------
    // Se la nuova configurazione è più probabile,
    // la mossa viene sempre accettata
    // --------------------------------------------------------

    if (log_ratio >= 0.0) {

        position = trial;

        return true;
    }


    // --------------------------------------------------------
    // Se la nuova configurazione è meno probabile,
    // viene accettata con probabilità exp(log_ratio)
    // --------------------------------------------------------

    const double acceptance =
        exp(log_ratio);


    if (rnd.Rannyu() <= acceptance) {

        position = trial;

        return true;
    }


    // Mossa rifiutata:
    // la configurazione rimane invariata
    return false;
}


// ============================================================
// Nome leggibile della transition probability
// ============================================================

string ProposalName(const string& proposal) {

    if (proposal == "uniform") {

        return "uniforme";
    }

    if (proposal == "gauss") {

        return "gaussiana";
    }

    return proposal;
}


// ============================================================
// Taratura automatica del passo della transition probability
//
// Per la proposal uniforme il parametro è Delta.
//
// Per la proposal gaussiana il parametro è sigma.
//
// Il parametro viene modificato attraverso brevi simulazioni
// preliminari fino a ottenere un'accettanza prossima al 50%.
//
// I dati generati durante questa fase non vengono utilizzati
// nelle misure finali.
// ============================================================

double TuneStep(const string& orbital,
                const string& proposal,
                Position start,
                double step,
                int steps_per_test,
                int max_iterations,
                Random& rnd) {

    const double target = 0.50;

    const double tolerance = 0.01;


    Position position = start;


    cout << endl;

    cout << "Taratura della proposal "
         << ProposalName(proposal)
         << " per l'orbitale "
         << orbital
         << endl;


    // ========================================================
    // Ciclo di tuning
    // ========================================================

    for (int iteration = 0;
         iteration < max_iterations;
         ++iteration) {

        int accepted = 0;


        // ----------------------------------------------------
        // Breve simulazione per stimare l'accettanza
        // ----------------------------------------------------

        for (int i = 0;
             i < steps_per_test;
             ++i) {

            if (MetropolisStep(
                    position,
                    step,
                    orbital,
                    proposal,
                    rnd)) {

                accepted++;
            }
        }


        const double acceptance =
            static_cast<double>(accepted)
            / steps_per_test;


        cout << "Iterazione "
             << iteration + 1
             << ": step = "
             << step
             << "   accettanza = "
             << acceptance
             << endl;


        // ----------------------------------------------------
        // Controllo del raggiungimento dell'accettanza target
        // ----------------------------------------------------

        if (fabs(acceptance - target) < tolerance) {

            cout << "Taratura completata." << endl;

            return step;
        }


        // ----------------------------------------------------
        // Accettanza troppo alta:
        // i passi sono troppo piccoli
        // ----------------------------------------------------

        if (acceptance > target) {

            step *= 1.10;
        }


        // ----------------------------------------------------
        // Accettanza troppo bassa:
        // i passi sono troppo grandi
        // ----------------------------------------------------

        else {

            step *= 0.90;
        }
    }


    cout << "Raggiunto il numero massimo di iterazioni "
         << "per la taratura."
         << endl;


    return step;
}


// ============================================================
// Studio dell'equilibrazione
//
// Questa funzione viene utilizzata solamente per la proposal
// uniforme, per produrre i grafici già inseriti nel notebook.
//
// I dati generati non vengono utilizzati nelle misure.
// ============================================================

void StudyEquilibration(const string& orbital,
                        const string& proposal,
                        Position start,
                        double step,
                        int n_steps,
                        const string& filename,
                        Random& rnd) {

    Position position = start;


    ofstream output(filename);


    if (!output.is_open()) {

        cerr << "Errore nell'apertura del file "
             << filename
             << endl;

        exit(1);
    }


    // --------------------------------------------------------
    // Configurazione iniziale
    // --------------------------------------------------------

    output
        << 0 << " "
        << Radius(position)
        << endl;


    // ========================================================
    // Evoluzione della catena
    // ========================================================

    for (int i = 1;
         i <= n_steps;
         ++i) {

        MetropolisStep(
            position,
            step,
            orbital,
            proposal,
            rnd
        );


        output
            << i << " "
            << Radius(position)
            << endl;
    }


    output.close();
}


// ============================================================
// Calcolo dell'incertezza statistica mediante data blocking
// ============================================================

double Error(double average,
             double average2,
             int n) {

    // Con un solo blocco non è possibile stimare l'errore
    if (n == 0) {

        return 0.0;
    }


    double variance =
        average2 - average * average;


    // Protezione da eventuali piccoli errori numerici
    if (variance < 0.0) {

        variance = 0.0;
    }


    return sqrt(variance / n);
}


// ============================================================
// Simulazione di produzione di un orbitale
//
// Per la proposal uniforme vengono mantenuti i nomi originali:
//
// radius_1s.dat
// positions_1s.dat
// series_1s.dat
//
// Per la proposal gaussiana viene aggiunto "_gauss":
//
// radius_1s_gauss.dat
// positions_1s_gauss.dat
// series_1s_gauss.dat
//
// In questo modo tutte le celle precedenti del notebook
// continuano a funzionare senza modifiche.
// ============================================================

SimulationResult SimulateOrbital(const string& orbital,
                                 const string& proposal,
                                 Position start,
                                 double step,
                                 int n_equil,
                                 int n_blocks,
                                 int steps_per_block,
                                 Random& rnd) {


    // ========================================================
    // Costruzione dei nomi dei file
    // ========================================================

    string suffix = "";


    // La simulazione uniforme conserva i nomi già utilizzati
    // nel notebook
    if (proposal == "gauss") {

        suffix = "_gauss";
    }


    const string radius_file =
        "radius_" + orbital + suffix + ".dat";

    const string positions_file =
        "positions_" + orbital + suffix + ".dat";

    const string series_file =
        "series_" + orbital + suffix + ".dat";


    ofstream radius_out(radius_file);

    ofstream positions_out(positions_file);

    ofstream series_out(series_file);


    if (!radius_out.is_open() ||
        !positions_out.is_open() ||
        !series_out.is_open()) {

        cerr << "Errore nell'apertura dei file di output."
             << endl;

        exit(1);
    }


    // ========================================================
    // Reset alla posizione iniziale
    // ========================================================

    Position position = start;


    // ========================================================
    // Equilibrazione
    //
    // I primi n_equil passi vengono scartati.
    // ========================================================

    int accepted_equil = 0;


    for (int i = 0;
         i < n_equil;
         ++i) {

        if (MetropolisStep(
                position,
                step,
                orbital,
                proposal,
                rnd)) {

            accepted_equil++;
        }
    }


    cout << endl;

    cout << "Orbitale: "
         << orbital
         << "   proposal: "
         << ProposalName(proposal)
         << endl;


    cout << "Accettanza durante l'equilibrazione = "
         << static_cast<double>(accepted_equil)
            / n_equil
         << endl;


    // ========================================================
    // Variabili per il data blocking
    // ========================================================

    double sum_prog = 0.0;

    double sum2_prog = 0.0;


    double final_average = 0.0;

    double final_error = 0.0;


    long long accepted = 0;

    long long attempted = 0;


    // Per il grafico 3D viene salvata una posizione
    // ogni 100 passi Metropolis
    const int position_stride = 100;


    long long global_step = 0;


    // ========================================================
    // Ciclo sui blocchi
    // ========================================================

    for (int block = 0;
         block < n_blocks;
         ++block) {

        double block_sum = 0.0;


        // ====================================================
        // Passi Metropolis all'interno del blocco
        // ====================================================

        for (int i = 0;
             i < steps_per_block;
             ++i) {


            // ------------------------------------------------
            // Evoluzione della catena
            // ------------------------------------------------

            if (MetropolisStep(
                    position,
                    step,
                    orbital,
                    proposal,
                    rnd)) {

                accepted++;
            }


            attempted++;


            // ------------------------------------------------
            // Misura della distanza dall'origine
            // ------------------------------------------------

            const double r =
                Radius(position);


            block_sum += r;


            // ------------------------------------------------
            // Serie completa di r
            // ------------------------------------------------

            series_out
                << global_step + 1 << " "
                << r
                << endl;


            // ------------------------------------------------
            // Sottoinsieme delle posizioni per il grafico 3D
            // ------------------------------------------------

            if (global_step % position_stride == 0) {

                positions_out
                    << position.x << " "
                    << position.y << " "
                    << position.z
                    << endl;
            }


            global_step++;
        }


        // ====================================================
        // Media del blocco
        // ====================================================

        const double block_average =
            block_sum / steps_per_block;


        // ====================================================
        // Somme progressive
        // ====================================================

        sum_prog +=
            block_average;


        sum2_prog +=
            block_average * block_average;


        // ====================================================
        // Media progressiva
        // ====================================================

        const double progressive_average =
            sum_prog / (block + 1);


        const double progressive_average2 =
            sum2_prog / (block + 1);


        // ====================================================
        // Errore progressivo
        // ====================================================

        const double error =
            Error(
                progressive_average,
                progressive_average2,
                block
            );


        // Salvataggio dei valori finali
        final_average =
            progressive_average;

        final_error =
            error;


        // ====================================================
        // File del data blocking
        //
        // Colonne:
        //
        // 1) numero del blocco
        // 2) numero cumulativo di passi
        // 3) media progressiva di <r>
        // 4) errore
        // ====================================================

        radius_out
            << block + 1 << " "
            << static_cast<long long>(block + 1)
               * steps_per_block << " "
            << progressive_average << " "
            << error
            << endl;
    }


    // ========================================================
    // Accettanza finale
    // ========================================================

    const double acceptance_rate =
        static_cast<double>(accepted)
        / attempted;


    cout << "Accettanza durante le misure = "
         << acceptance_rate
         << endl;


    cout << "Valore finale di <r> = "
         << final_average
         << " +/- "
         << final_error
         << endl;


    // ========================================================
    // Chiusura dei file
    // ========================================================

    radius_out.close();

    positions_out.close();

    series_out.close();


    // ========================================================
    // Restituzione dei risultati finali
    // ========================================================

    SimulationResult result;

    result.mean =
        final_average;

    result.error =
        final_error;

    result.acceptance =
        acceptance_rate;


    return result;
}


// ============================================================
// MAIN
// ============================================================

int main() {


    // ========================================================
    // Inizializzazione del generatore di numeri casuali
    // ========================================================

    Random rnd;


    int seed[4];

    int p1;
    int p2;


    // --------------------------------------------------------
    // Lettura dei numeri primi
    // --------------------------------------------------------

    ifstream primes("Primes");


    if (primes.is_open()) {

        primes >> p1 >> p2;
    }

    else {

        cerr << "PROBLEMA: impossibile aprire il file Primes"
             << endl;

        return 1;
    }


    primes.close();


    // --------------------------------------------------------
    // Lettura del seme iniziale
    // --------------------------------------------------------

    ifstream input("seed.in");


    string property;


    if (input.is_open()) {

        while (input >> property) {

            if (property == "RANDOMSEED") {

                input
                    >> seed[0]
                    >> seed[1]
                    >> seed[2]
                    >> seed[3];


                rnd.SetRandom(
                    seed,
                    p1,
                    p2
                );
            }
        }
    }

    else {

        cerr << "PROBLEMA: impossibile aprire il file seed.in"
             << endl;

        return 1;
    }


    input.close();


    // ========================================================
    // Parametri generali
    // ========================================================

    const int n_blocks =
        100;


    const int steps_per_block =
        10000;


    // Numero totale di misure:
    //
    // M = 100 * 10000 = 10^6


    // Passi iniziali scartati
    const int n_equil =
        10000;


    // ========================================================
    // Parametri del tuning
    // ========================================================

    const int tuning_steps =
        5000;


    const int tuning_iterations =
        30;


    // ========================================================
    // Passi utilizzati solamente per i grafici
    // dell'equilibrazione
    // ========================================================

    const int equilibration_study_steps =
        30000;


    // ========================================================
    // Posizioni iniziali
    // ========================================================

    Position start_1s =
        {0.0, 0.0, 1.0};


    // Per il 2p scegliamo una posizione fuori
    // dal piano nodale z = 0
    Position start_2p =
        {0.0, 0.0, 2.0};


    // ========================================================
    // ========================================================
    //
    // PROPOSAL UNIFORME
    //
    // ========================================================
    // ========================================================


    // ========================================================
    // Stato 1s - uniforme
    // ========================================================

    double delta_1s =
        1.0;


    delta_1s =
        TuneStep(
            "1s",
            "uniform",
            start_1s,
            delta_1s,
            tuning_steps,
            tuning_iterations,
            rnd
        );


    cout << "Delta finale uniforme per 1s = "
         << delta_1s
         << endl;


    // ========================================================
    // Studio dell'equilibrazione dello stato 1s
    //
    // Manteniamo solamente la proposal uniforme perché
    // questi grafici sono già utilizzati nel notebook.
    // ========================================================

    Position equil_start_near_1s =
        {0.0, 0.0, 1.0};


    StudyEquilibration(
        "1s",
        "uniform",
        equil_start_near_1s,
        delta_1s,
        equilibration_study_steps,
        "equilibration_1s_near.dat",
        rnd
    );


    Position equil_start_far_1s =
        {1000.0, 0.0, 0.0};


    StudyEquilibration(
        "1s",
        "uniform",
        equil_start_far_1s,
        delta_1s,
        equilibration_study_steps,
        "equilibration_1s_far.dat",
        rnd
    );


    const SimulationResult uniform_1s =
        SimulateOrbital(
            "1s",
            "uniform",
            start_1s,
            delta_1s,
            n_equil,
            n_blocks,
            steps_per_block,
            rnd
        );


    // ========================================================
    // Stato 2p - uniforme
    // ========================================================

    double delta_2p =
        2.0;


    delta_2p =
        TuneStep(
            "2p",
            "uniform",
            start_2p,
            delta_2p,
            tuning_steps,
            tuning_iterations,
            rnd
        );


    cout << "Delta finale uniforme per 2p = "
         << delta_2p
         << endl;


    const SimulationResult uniform_2p =
        SimulateOrbital(
            "2p",
            "uniform",
            start_2p,
            delta_2p,
            n_equil,
            n_blocks,
            steps_per_block,
            rnd
        );


    // ========================================================
    // ========================================================
    //
    // PROPOSAL GAUSSIANA
    //
    // ========================================================
    // ========================================================


    // ========================================================
    // Stato 1s - gaussiana
    // ========================================================

    // Valore iniziale di prova per sigma
    double sigma_1s =
        1.0;


    sigma_1s =
        TuneStep(
            "1s",
            "gauss",
            start_1s,
            sigma_1s,
            tuning_steps,
            tuning_iterations,
            rnd
        );


    cout << "Sigma finale gaussiana per 1s = "
         << sigma_1s
         << endl;


    const SimulationResult gauss_1s =
        SimulateOrbital(
            "1s",
            "gauss",
            start_1s,
            sigma_1s,
            n_equil,
            n_blocks,
            steps_per_block,
            rnd
        );


    // ========================================================
    // Stato 2p - gaussiana
    // ========================================================

    double sigma_2p =
        2.0;


    sigma_2p =
        TuneStep(
            "2p",
            "gauss",
            start_2p,
            sigma_2p,
            tuning_steps,
            tuning_iterations,
            rnd
        );


    cout << "Sigma finale gaussiana per 2p = "
         << sigma_2p
         << endl;


    const SimulationResult gauss_2p =
        SimulateOrbital(
            "2p",
            "gauss",
            start_2p,
            sigma_2p,
            n_equil,
            n_blocks,
            steps_per_block,
            rnd
        );


    // ========================================================
    // Confronto finale tra proposal uniforme e gaussiana
    // ========================================================

    cout << endl;
    cout << "============================================="
         << endl;

    cout << "CONFRONTO FINALE"
         << endl;

    cout << "============================================="
         << endl;


    cout << endl;

    cout << "Stato 1s"
         << endl;


    cout << "Uniforme : <r> = "
         << uniform_1s.mean
         << " +/- "
         << uniform_1s.error
         << "   accettanza = "
         << uniform_1s.acceptance
         << endl;


    cout << "Gaussiana: <r> = "
         << gauss_1s.mean
         << " +/- "
         << gauss_1s.error
         << "   accettanza = "
         << gauss_1s.acceptance
         << endl;


    cout << "Valore teorico = 1.5"
         << endl;


    cout << endl;

    cout << "Stato 2p"
         << endl;


    cout << "Uniforme : <r> = "
         << uniform_2p.mean
         << " +/- "
         << uniform_2p.error
         << "   accettanza = "
         << uniform_2p.acceptance
         << endl;


    cout << "Gaussiana: <r> = "
         << gauss_2p.mean
         << " +/- "
         << gauss_2p.error
         << "   accettanza = "
         << gauss_2p.acceptance
         << endl;


    cout << "Valore teorico = 5.0"
         << endl;


    cout << endl;
    cout << "============================================="
         << endl;


    // ========================================================
    // Salvataggio del seme finale
    // ========================================================

    rnd.SaveSeed();


    return 0;
}