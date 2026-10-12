#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <string>

#include "random.h"

using namespace std;


// ============================================================
// Struttura utilizzata per restituire i risultati
// di una simulazione Variational Monte Carlo.
// ============================================================

struct VMCResult {

    double energy;
    double error;
    double acceptance;
};


// ============================================================
// Funzione d'onda di prova
//
// Psi_T(x) = exp[-(x-mu)^2 / (2 sigma^2)]
//          + exp[-(x+mu)^2 / (2 sigma^2)]
//
// La normalizzazione non è necessaria perché nel Metropolis
// vengono utilizzati solamente rapporti di probabilità.
// ============================================================

double psi(double x, double mu, double sigma) {

    double sigma2 = sigma * sigma;

    double g_minus =
        exp(-(x - mu) * (x - mu) / (2.0 * sigma2));

    double g_plus =
        exp(-(x + mu) * (x + mu) / (2.0 * sigma2));

    return g_minus + g_plus;
}


// ============================================================
// Distribuzione di probabilità campionata dal Metropolis
//
// p(x) proporzionale a |Psi_T(x)|^2
// ============================================================

double probability(double x, double mu, double sigma) {

    double psi_value = psi(x, mu, sigma);

    return psi_value * psi_value;
}


// ============================================================
// Potenziale esterno
//
// V(x) = x^4 - (5/2)x^2
// ============================================================

double potential(double x) {

    return pow(x, 4) - 2.5 * x * x;
}


// ============================================================
// Energia locale
//
// E_loc(x) = H Psi_T(x) / Psi_T(x)
//
// Con hbar = 1 e m = 1:
//
// E_loc(x) = -1/2 Psi_T''(x)/Psi_T(x) + V(x)
// ============================================================

double local_energy(double x, double mu, double sigma) {

    double sigma2 = sigma * sigma;
    double sigma4 = sigma2 * sigma2;


    // Prima gaussiana centrata in +mu
    double g_minus =
        exp(-(x - mu) * (x - mu) / (2.0 * sigma2));


    // Seconda gaussiana centrata in -mu
    double g_plus =
        exp(-(x + mu) * (x + mu) / (2.0 * sigma2));


    // Derivata seconda della funzione d'onda di prova
    double psi_second =
        (
            ((x - mu) * (x - mu) / sigma4 - 1.0 / sigma2)
            * g_minus
        )
        +
        (
            ((x + mu) * (x + mu) / sigma4 - 1.0 / sigma2)
            * g_plus
        );


    double psi_value =
        g_minus + g_plus;


    // Energia cinetica locale
    double kinetic =
        -0.5 * psi_second / psi_value;


    // Energia locale totale
    return kinetic + potential(x);
}


// ============================================================
// Singolo passo Metropolis nello spazio delle posizioni.
//
// Viene proposta:
//
// x_new = x_old + U(-delta_x, delta_x)
//
// Essendo la probabilità di transizione simmetrica:
//
// A = min[1, p(x_new)/p(x_old)]
// ============================================================

bool metropolis_step(Random &rnd,
                     double &x,
                     double mu,
                     double sigma,
                     double delta_x) {

    // Proposta della nuova posizione
    double x_new =
        x + rnd.Rannyu(-delta_x, delta_x);


    // Probabilità della configurazione vecchia e nuova
    double p_old =
        probability(x, mu, sigma);

    double p_new =
        probability(x_new, mu, sigma);


    // Probabilità di accettazione
    double acceptance =
        min(1.0, p_new / p_old);


    // Accettazione o rifiuto della nuova posizione
    if (rnd.Rannyu() < acceptance) {

        x = x_new;

        return true;
    }


    return false;
}


// ============================================================
// Errore statistico sulla media progressiva.
//
// mean  = <A>
// mean2 = <A^2>
// ============================================================

double error(double mean,
             double mean2,
             int n) {

    if (n <= 1)
        return 0.0;


    double variance =
        mean2 - mean * mean;


    // Protezione da piccoli valori negativi dovuti
    // agli errori numerici di arrotondamento.
    variance = max(0.0, variance);


    return sqrt(
        variance
        / static_cast<double>(n - 1)
    );
}


// ============================================================
// Calcolo dell'energia variazionale mediante VMC.
//
// Questa funzione viene richiamata dal Simulated Annealing
// ogni volta che deve essere valutata una nuova coppia
// di parametri (mu, sigma).
//
// Restituisce:
//  - energia media
//  - errore statistico
//  - accettazione del Metropolis
// ============================================================

VMCResult vmc_energy(Random &rnd,
                     double mu,
                     double sigma,
                     double delta_x,
                     int n_blocks,
                     int steps_per_block,
                     int equilibration_steps) {

    // Posizione iniziale della catena
    double x = 0.0;


    // --------------------------------------------------------
    // Equilibrazione
    // --------------------------------------------------------

    for (int i = 0; i < equilibration_steps; ++i) {

        metropolis_step(
            rnd,
            x,
            mu,
            sigma,
            delta_x
        );
    }


    double progressive_sum = 0.0;
    double progressive_sum2 = 0.0;


    long long accepted = 0;
    long long attempted = 0;


    // --------------------------------------------------------
    // Data blocking
    // --------------------------------------------------------

    for (int block = 0; block < n_blocks; ++block) {

        double block_sum = 0.0;


        for (int i = 0; i < steps_per_block; ++i) {


            bool move_accepted =
                metropolis_step(
                    rnd,
                    x,
                    mu,
                    sigma,
                    delta_x
                );


            if (move_accepted)
                accepted++;


            attempted++;


            // Misura dell'energia locale
            block_sum +=
                local_energy(
                    x,
                    mu,
                    sigma
                );
        }


        // Media dell'energia nel blocco corrente
        double block_average =
            block_sum
            / static_cast<double>(steps_per_block);


        progressive_sum +=
            block_average;


        progressive_sum2 +=
            block_average * block_average;
    }


    // --------------------------------------------------------
    // Energia media finale
    // --------------------------------------------------------

    double mean =
        progressive_sum
        / static_cast<double>(n_blocks);


    double mean2 =
        progressive_sum2
        / static_cast<double>(n_blocks);


    double uncertainty =
        error(
            mean,
            mean2,
            n_blocks
        );


    double acceptance_rate =
        static_cast<double>(accepted)
        / static_cast<double>(attempted);


    return {
        mean,
        uncertainty,
        acceptance_rate
    };
}


// ============================================================
// Simulazione VMC finale ad alta statistica.
//
// Utilizza i parametri ottimizzati dal Simulated Annealing
// e salva:
//
//  - energia progressiva
//  - errore progressivo
//  - posizioni x campionate
// ============================================================

VMCResult final_vmc(Random &rnd,
                    double mu,
                    double sigma,
                    double delta_x,
                    int n_blocks,
                    int steps_per_block,
                    int equilibration_steps,
                    int sample_stride) {


    // --------------------------------------------------------
    // File contenente l'energia progressiva
    // --------------------------------------------------------

    ofstream energy_file("final_energy.txt");


    energy_file
        << "# Energia variazionale finale ottenuta tramite VMC\n"
        << "# I parametri mu e sigma sono quelli ottimizzati dal Simulated Annealing\n"
        << "# mu = " << mu << "\n"
        << "# sigma = " << sigma << "\n"
        << "# delta_x = " << delta_x << "\n"
        << "# Numero blocchi = " << n_blocks << "\n"
        << "# Passi per blocco = " << steps_per_block << "\n"
        << "# Colonne:\n"
        << "# 1 = numero del blocco\n"
        << "# 2 = energia media del singolo blocco\n"
        << "# 3 = energia media progressiva\n"
        << "# 4 = incertezza statistica progressiva\n"
        << "# block block_energy progressive_energy error\n";


    // --------------------------------------------------------
    // File contenente le posizioni campionate
    // --------------------------------------------------------

    ofstream sample_file("final_samples.txt");


    sample_file
        << "# Posizioni x campionate con il Variational Monte Carlo finale\n"
        << "# La distribuzione campionata e' proporzionale a |Psi_T(x)|^2\n"
        << "# mu = " << mu << "\n"
        << "# sigma = " << sigma << "\n"
        << "# delta_x = " << delta_x << "\n"
        << "# Viene salvato un punto ogni "
        << sample_stride
        << " passi Metropolis\n"
        << "# Colonne:\n"
        << "# 1 = posizione x campionata\n"
        << "# x\n";


    // Posizione iniziale
    double x = 0.0;


    // --------------------------------------------------------
    // Equilibrazione
    // --------------------------------------------------------

    for (int i = 0; i < equilibration_steps; ++i) {

        metropolis_step(
            rnd,
            x,
            mu,
            sigma,
            delta_x
        );
    }


    double progressive_sum = 0.0;
    double progressive_sum2 = 0.0;


    long long accepted = 0;
    long long attempted = 0;

    long long total_step = 0;


    // --------------------------------------------------------
    // Ciclo sui blocchi
    // --------------------------------------------------------

    for (int block = 0; block < n_blocks; ++block) {

        double block_sum = 0.0;


        for (int i = 0; i < steps_per_block; ++i) {


            bool move_accepted =
                metropolis_step(
                    rnd,
                    x,
                    mu,
                    sigma,
                    delta_x
                );


            if (move_accepted)
                accepted++;


            attempted++;

            total_step++;


            // Misura dell'energia locale
            block_sum +=
                local_energy(
                    x,
                    mu,
                    sigma
                );


            // Salvataggio periodico delle configurazioni
            if (total_step % sample_stride == 0) {

                sample_file << x << '\n';
            }
        }


        // ----------------------------------------------------
        // Media del blocco
        // ----------------------------------------------------

        double block_average =
            block_sum
            / static_cast<double>(steps_per_block);


        progressive_sum +=
            block_average;


        progressive_sum2 +=
            block_average * block_average;


        int n = block + 1;


        // ----------------------------------------------------
        // Media progressiva
        // ----------------------------------------------------

        double progressive_average =
            progressive_sum
            / static_cast<double>(n);


        double progressive_average2 =
            progressive_sum2
            / static_cast<double>(n);


        double uncertainty =
            error(
                progressive_average,
                progressive_average2,
                n
            );


        // Salvataggio sul file
        energy_file
            << n << " "
            << block_average << " "
            << progressive_average << " "
            << uncertainty
            << '\n';
    }


    // --------------------------------------------------------
    // Risultato finale
    // --------------------------------------------------------

    double final_energy =
        progressive_sum
        / static_cast<double>(n_blocks);


    double final_energy2 =
        progressive_sum2
        / static_cast<double>(n_blocks);


    double final_error =
        error(
            final_energy,
            final_energy2,
            n_blocks
        );


    double acceptance_rate =
        static_cast<double>(accepted)
        / static_cast<double>(attempted);


    energy_file.close();

    sample_file.close();


    return {
        final_energy,
        final_error,
        acceptance_rate
    };
}


// ============================================================
// MAIN
// ============================================================

int main() {


    // ========================================================
    // INIZIALIZZAZIONE DEL GENERATORE DI NUMERI CASUALI
    // ========================================================

    Random rnd;


    int seed[4];

    int p1, p2;


    // Lettura dei numeri primi
    ifstream primes("Primes");


    if (primes.is_open()) {

        primes >> p1 >> p2;

        primes.close();
    }

    else {

        cerr
            << "Errore: impossibile aprire il file Primes."
            << endl;

        return 1;
    }


    // Lettura del seed
    ifstream input("seed.in");

    string property;


    if (input.is_open()) {

        while (!input.eof()) {

            input >> property;


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

        input.close();
    }

    else {

        cerr
            << "Errore: impossibile aprire il file seed.in."
            << endl;

        return 1;
    }


    // ========================================================
    // PARAMETRI DEL METROPOLIS
    // ========================================================

    // Durante il Simulated Annealing viene utilizzato il
    // valore testato inizialmente, che fornisce una buona
    // accettazione per i parametri di partenza.
    const double delta_x_sa = 1.2;


    // Nella simulazione finale viene utilizzata una proposta
    // più ampia perché la distribuzione ottimizzata risulta
    // più larga e con delta_x = 1.2 l'accettazione era circa 73%.
    const double delta_x_final = 2.5;


    // ========================================================
    // STATISTICA DEL VMC DURANTE IL SIMULATED ANNEALING
    // ========================================================

    // Ogni valutazione della funzione costo utilizza:
    //
    // 100 blocchi x 2000 passi = 200000 misure

    const int vmc_blocks_sa = 100;

    const int vmc_steps_sa = 2000;

    const int vmc_equilibration_sa = 5000;


    // ========================================================
    // PARAMETRI DEL SIMULATED ANNEALING
    // ========================================================

    // Parametri variazionali iniziali
    double mu = 1.0;

    double sigma = 0.5;


    // Temperatura fittizia iniziale
    double temperature = 0.10;


    // Raffreddamento geometrico:
    //
    // T_(n+1) = cooling * T_n
    const double cooling = 0.98;


    // Numero di passi del Simulated Annealing
    const int n_sa_steps = 200;


    // Ampiezza delle proposte nello spazio dei parametri
    const double delta_mu = 0.10;

    const double delta_sigma = 0.10;


    // Sigma rappresenta una larghezza e deve essere positivo
    const double sigma_min = 0.05;


    // ========================================================
    // FILE DEL SIMULATED ANNEALING
    // ========================================================

    ofstream annealing_file("annealing.txt");


    annealing_file
        << "# Evoluzione del Simulated Annealing\n"
        << "# La funzione costo e' l'energia variazionale <H>\n"
        << "# Parametri iniziali:\n"
        << "# mu = " << mu << "\n"
        << "# sigma = " << sigma << "\n"
        << "# Temperatura iniziale = " << temperature << "\n"
        << "# Cooling factor = " << cooling << "\n"
        << "# delta_mu = " << delta_mu << "\n"
        << "# delta_sigma = " << delta_sigma << "\n"
        << "# delta_x utilizzato nel VMC = " << delta_x_sa << "\n"
        << "# Colonne:\n"
        << "# 1 = passo del Simulated Annealing\n"
        << "# 2 = temperatura T\n"
        << "# 3 = parametro mu della configurazione corrente\n"
        << "# 4 = parametro sigma della configurazione corrente\n"
        << "# 5 = energia variazionale della configurazione corrente\n"
        << "# 6 = errore statistico sull'energia\n"
        << "# 7 = mossa SA accettata: 1 si, 0 no\n"
        << "# 8 = migliore energia incontrata fino a quel passo\n"
        << "# step T mu sigma energy error accepted best_energy\n";


    // ========================================================
    // ENERGIA DELLA CONFIGURAZIONE INIZIALE
    // ========================================================

    cout << endl;

    cout
        << "Calcolo energia iniziale..."
        << endl;


    VMCResult current =
        vmc_energy(
            rnd,
            mu,
            sigma,
            delta_x_sa,
            vmc_blocks_sa,
            vmc_steps_sa,
            vmc_equilibration_sa
        );


    // --------------------------------------------------------
    // Migliore configurazione incontrata
    // --------------------------------------------------------

    double best_mu =
        mu;


    double best_sigma =
        sigma;


    double best_energy =
        current.energy;


    double best_error =
        current.error;


    // Salvataggio della configurazione iniziale
    annealing_file
        << 0 << " "
        << temperature << " "
        << mu << " "
        << sigma << " "
        << current.energy << " "
        << current.error << " "
        << 1 << " "
        << best_energy
        << '\n';


    cout
        << "Energia iniziale = "
        << current.energy
        << " +/- "
        << current.error
        << endl;


    // Contatore delle mosse SA accettate
    int sa_accepted = 0;


    // ========================================================
    // CICLO DI SIMULATED ANNEALING
    // ========================================================

    for (int step = 1;
         step <= n_sa_steps;
         ++step) {


        // ----------------------------------------------------
        // Proposta di una nuova coppia (mu, sigma)
        // ----------------------------------------------------

        double mu_new =
            mu
            + rnd.Rannyu(
                -delta_mu,
                delta_mu
            );


        double sigma_new =
            sigma
            + rnd.Rannyu(
                -delta_sigma,
                delta_sigma
            );


        bool accepted = false;


        // ----------------------------------------------------
        // Sigma deve rimanere positivo
        // ----------------------------------------------------

        if (sigma_new > sigma_min) {


            // ------------------------------------------------
            // Calcolo dell'energia della nuova configurazione
            // ------------------------------------------------

            VMCResult candidate =
                vmc_energy(
                    rnd,
                    mu_new,
                    sigma_new,
                    delta_x_sa,
                    vmc_blocks_sa,
                    vmc_steps_sa,
                    vmc_equilibration_sa
                );


            // Differenza di energia
            double delta_energy =
                candidate.energy
                - current.energy;


            // ------------------------------------------------
            // Regola di accettazione del Simulated Annealing
            //
            // DeltaE <= 0:
            //      accettazione certa
            //
            // DeltaE > 0:
            //
            //      P = exp(-DeltaE / T)
            // ------------------------------------------------

            if (delta_energy <= 0.0) {

                accepted = true;
            }

            else {

                double sa_probability =
                    exp(
                        -delta_energy
                        / temperature
                    );


                if (rnd.Rannyu() < sa_probability) {

                    accepted = true;
                }
            }


            // ------------------------------------------------
            // Aggiornamento della configurazione corrente
            // ------------------------------------------------

            if (accepted) {

                mu =
                    mu_new;


                sigma =
                    sigma_new;


                current =
                    candidate;


                sa_accepted++;


                // --------------------------------------------
                // Aggiornamento della migliore configurazione
                // incontrata fino a questo punto
                // --------------------------------------------

                if (current.energy < best_energy) {

                    best_energy =
                        current.energy;


                    best_error =
                        current.error;


                    best_mu =
                        mu;


                    best_sigma =
                        sigma;
                }
            }
        }


        // ----------------------------------------------------
        // Salvataggio dello stato corrente
        // ----------------------------------------------------

        annealing_file
            << step << " "
            << temperature << " "
            << mu << " "
            << sigma << " "
            << current.energy << " "
            << current.error << " "
            << accepted << " "
            << best_energy
            << '\n';


        // ----------------------------------------------------
        // Raffreddamento
        //
        // T_(n+1) = cooling * T_n
        // ----------------------------------------------------

        temperature *=
            cooling;


        // ----------------------------------------------------
        // Informazioni periodiche a terminale
        // ----------------------------------------------------

        if (step % 20 == 0) {

            cout
                << "Step "
                << step
                << " / "
                << n_sa_steps
                << "   T = "
                << temperature
                << "   E = "
                << current.energy
                << "   mu = "
                << mu
                << "   sigma = "
                << sigma
                << endl;
        }
    }


    annealing_file.close();


    // ========================================================
    // RISULTATO DEL SIMULATED ANNEALING
    // ========================================================

    double sa_acceptance =
        static_cast<double>(sa_accepted)
        / static_cast<double>(n_sa_steps);


    cout << endl;


    cout
        << "Simulated Annealing terminato"
        << endl;


    cout
        << "----------------------------"
        << endl;


    cout
        << "Migliore mu    = "
        << best_mu
        << endl;


    cout
        << "Migliore sigma = "
        << best_sigma
        << endl;


    cout
        << "Energia SA     = "
        << best_energy
        << " +/- "
        << best_error
        << endl;


    cout
        << "Accettazione SA = "
        << 100.0 * sa_acceptance
        << " %"
        << endl;


    // ========================================================
    // SIMULAZIONE FINALE AD ALTA STATISTICA
    // ========================================================

    // I parametri migliori trovati dal Simulated Annealing
    // vengono utilizzati per una nuova simulazione VMC con
    // statistica maggiore.

    const int final_blocks =
        1000;


    const int final_steps_per_block =
        1000;


    const int final_equilibration =
        10000;


    // Viene salvata una configurazione ogni 10 passi
    const int sample_stride =
        10;


    cout << endl;


    cout
        << "Simulazione finale con i parametri ottimizzati..."
        << endl;


    VMCResult final_result =
        final_vmc(
            rnd,
            best_mu,
            best_sigma,
            delta_x_final,
            final_blocks,
            final_steps_per_block,
            final_equilibration,
            sample_stride
        );


    // ========================================================
    // RISULTATO FINALE
    // ========================================================

    cout << endl;


    cout
        << "Risultato finale"
        << endl;


    cout
        << "---------------"
        << endl;


    cout
        << "mu       = "
        << best_mu
        << endl;


    cout
        << "sigma    = "
        << best_sigma
        << endl;


    cout
        << "<H>      = "
        << final_result.energy
        << " +/- "
        << final_result.error
        << endl;


    cout
        << "Accettazione Metropolis = "
        << 100.0 * final_result.acceptance
        << " %"
        << endl;


    // ========================================================
    // FILE RIASSUNTIVO DEI PARAMETRI OTTIMIZZATI
    // ========================================================

    ofstream best_file(
        "best_parameters.txt"
    );


    best_file
        << "# Risultato dell'ottimizzazione variazionale\n"
        << "# I parametri sono stati ottenuti tramite Simulated Annealing\n"
        << "# L'energia finale e' stata ricalcolata con VMC ad alta statistica\n"
        << "# Numero blocchi finali = " << final_blocks << "\n"
        << "# Passi per blocco finali = " << final_steps_per_block << "\n"
        << "# delta_x finale = " << delta_x_final << "\n"
        << "# Colonne:\n"
        << "# 1 = mu ottimizzato\n"
        << "# 2 = sigma ottimizzato\n"
        << "# 3 = energia variazionale finale\n"
        << "# 4 = errore statistico finale\n"
        << "# 5 = accettazione Metropolis finale\n"
        << "# mu sigma energy error acceptance\n"
        << best_mu << " "
        << best_sigma << " "
        << final_result.energy << " "
        << final_result.error << " "
        << final_result.acceptance
        << '\n';


    best_file.close();


    // ========================================================
    // SALVATAGGIO DEL SEED
    // ========================================================

    rnd.SaveSeed();


    return 0;
}