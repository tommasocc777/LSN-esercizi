#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <string>

#include "random.h"

using namespace std;


// ============================================================
// Funzione d'onda di prova
//
// Psi_T(x) = exp[-(x-mu)^2 / (2 sigma^2)]
//          + exp[-(x+mu)^2 / (2 sigma^2)]
//
// La normalizzazione non è necessaria perché nel Metropolis
// utilizziamo solamente rapporti di probabilità.
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
// Distribuzione di probabilità da campionare con Metropolis
//
// p(x) proporzionale a |Psi_T(x)|^2.
//
// Non serve conoscere la costante di normalizzazione perché
// si semplifica nel rapporto p(x_new)/p(x_old).
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
//
// La derivata seconda delle due gaussiane viene calcolata
// analiticamente.
// ============================================================

double local_energy(double x, double mu, double sigma) {

    double sigma2 = sigma * sigma;
    double sigma4 = sigma2 * sigma2;


    // Prima gaussiana, centrata in +mu
    double g_minus =
        exp(-(x - mu) * (x - mu) / (2.0 * sigma2));


    // Seconda gaussiana, centrata in -mu
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


    double psi_value = g_minus + g_plus;


    // Contributo cinetico:
    // -1/2 * Psi''(x) / Psi(x)
    double kinetic =
        -0.5 * psi_second / psi_value;


    // Energia locale totale:
    // energia cinetica + potenziale
    return kinetic + potential(x);
}


// ============================================================
// Singolo passo dell'algoritmo di Metropolis.
//
// La nuova posizione viene proposta uniformemente:
//
// x_new = x_old + U(-delta, delta)
//
// Poiché la probabilità di transizione è simmetrica,
// l'accettazione è:
//
// A = min(1, p(x_new)/p(x_old))
//
// La funzione restituisce true se il passo viene accettato,
// false se viene rifiutato.
// ============================================================

bool metropolis_step(Random &rnd,
                     double &x,
                     double mu,
                     double sigma,
                     double delta) {

    // Proposta di una nuova posizione
    double x_new =
        x + rnd.Rannyu(-delta, delta);


    // Probabilità della vecchia e della nuova configurazione
    double p_old =
        probability(x, mu, sigma);

    double p_new =
        probability(x_new, mu, sigma);


    // Probabilità di accettazione Metropolis
    double acceptance =
        min(1.0, p_new / p_old);


    // Accetta o rifiuta il nuovo punto
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
//
// L'errore viene calcolato utilizzando le medie dei blocchi.
// ============================================================

double error(double mean,
             double mean2,
             int n) {

    if (n <= 1)
        return 0.0;

    return sqrt(
        (mean2 - mean * mean)
        / static_cast<double>(n - 1)
    );
}


// ============================================================
// MAIN
// ============================================================

int main(int argc, char *argv[]) {


    // --------------------------------------------------------
    // PARAMETRI VARIAZIONALI
    //
    // Per ora utilizziamo valori scelti manualmente.
    // Nell'esercizio 08.2 saranno ottimizzati automaticamente
    // tramite Simulated Annealing.
    //
    // È possibile passarli anche da terminale:
    //
    // ./main.exe mu sigma
    // --------------------------------------------------------

    double mu    = 0.8;
    double sigma = 0.6;


    if (argc >= 3) {

        mu    = stod(argv[1]);
        sigma = stod(argv[2]);
    }


    // Sigma deve essere strettamente positivo
    if (sigma <= 0.0) {

        cerr << "Errore: sigma deve essere positivo." << endl;

        return 1;
    }


    // --------------------------------------------------------
    // PARAMETRI DELLA SIMULAZIONE MONTE CARLO
    // --------------------------------------------------------

    const int n_blocks = 1000;

    const int steps_per_block = 1000;

    const int equilibration_steps = 2000;


    // Semilarghezza della distribuzione uniforme utilizzata
    // per proporre i nuovi punti nel Metropolis.
    //
    // Il valore va controllato tramite il tasso di accettazione.
    const double delta = 1.2;


    // Salviamo una posizione ogni sample_stride passi
    // per evitare un file eccessivamente grande.
    const int sample_stride = 10;


    // --------------------------------------------------------
    // INIZIALIZZAZIONE DEL GENERATORE DI NUMERI CASUALI
    // --------------------------------------------------------

    Random rnd;

    int seed[4];

    int p1, p2;


    ifstream primes("Primes");

    if (primes.is_open()) {

        primes >> p1 >> p2;

        primes.close();
    }
    else {

        cerr << "Errore: impossibile aprire il file Primes." << endl;

        return 1;
    }


    ifstream input("seed.in");

    string property;


    if (input.is_open()) {

        while (!input.eof()) {

            input >> property;


            if (property == "RANDOMSEED") {

                input >> seed[0]
                      >> seed[1]
                      >> seed[2]
                      >> seed[3];

                rnd.SetRandom(seed, p1, p2);
            }
        }

        input.close();
    }
    else {

        cerr << "Errore: impossibile aprire il file seed.in." << endl;

        return 1;
    }


    // --------------------------------------------------------
    // FILE DI OUTPUT PER L'ENERGIA
    // --------------------------------------------------------

    ofstream energy_file("energy.txt");


    energy_file
        << "# Energia variazionale calcolata con Variational Monte Carlo\n"
        << "# Parametri: mu = " << mu
        << ", sigma = " << sigma << "\n"
        << "# Ogni riga corrisponde a un blocco della simulazione\n"
        << "# Colonne:\n"
        << "# 1 = numero del blocco\n"
        << "# 2 = energia media del singolo blocco\n"
        << "# 3 = energia media progressiva <H>\n"
        << "# 4 = incertezza statistica sulla media progressiva\n"
        << "# block  block_energy  progressive_energy  error\n";


    // --------------------------------------------------------
    // FILE DI OUTPUT DELLE POSIZIONI CAMPIONATE
    // --------------------------------------------------------

    ofstream sample_file("samples.txt");


    sample_file
        << "# Posizioni x campionate con l'algoritmo di Metropolis\n"
        << "# La distribuzione campionata e' proporzionale a |Psi_T(x)|^2\n"
        << "# Viene salvato un punto ogni "
        << sample_stride
        << " passi Metropolis\n"
        << "# Parametri: mu = " << mu
        << ", sigma = " << sigma << "\n"
        << "# Colonne:\n"
        << "# 1 = posizione x campionata\n"
        << "# x\n";


    // --------------------------------------------------------
    // POSIZIONE INIZIALE DELLA CATENA DI MARKOV
    // --------------------------------------------------------

    double x = 0.0;


    // --------------------------------------------------------
    // EQUILIBRAZIONE
    //
    // Scartiamo i primi passi della catena perché la
    // distribuzione iniziale non è ancora rappresentativa
    // della distribuzione target |Psi_T(x)|^2.
    // --------------------------------------------------------

    for (int i = 0; i < equilibration_steps; ++i) {

        metropolis_step(
            rnd,
            x,
            mu,
            sigma,
            delta
        );
    }


    // --------------------------------------------------------
    // DATA BLOCKING
    // --------------------------------------------------------

    double progressive_sum  = 0.0;

    double progressive_sum2 = 0.0;


    // Contatori per il tasso di accettazione Metropolis
    long long accepted  = 0;

    long long attempted = 0;


    // Contatore globale dei passi della simulazione
    long long total_step = 0;


    // Ciclo sui blocchi
    for (int block = 0; block < n_blocks; ++block) {

        double block_sum = 0.0;


        // Ciclo sui passi del singolo blocco
        for (int i = 0; i < steps_per_block; ++i) {


            // Esegue un passo Metropolis
            bool move_accepted =
                metropolis_step(
                    rnd,
                    x,
                    mu,
                    sigma,
                    delta
                );


            if (move_accepted)
                accepted++;


            attempted++;

            total_step++;


            // Calcolo dell'energia locale nella posizione
            // attualmente visitata dalla catena
            block_sum +=
                local_energy(
                    x,
                    mu,
                    sigma
                );


            // Salvataggio periodico delle posizioni campionate
            if (total_step % sample_stride == 0) {

                sample_file << x << '\n';
            }
        }


        // ----------------------------------------------------
        // Media dell'energia nel blocco corrente
        // ----------------------------------------------------

        double block_average =
            block_sum
            / static_cast<double>(steps_per_block);


        // Aggiornamento delle somme progressive
        progressive_sum +=
            block_average;

        progressive_sum2 +=
            block_average * block_average;


        int n = block + 1;


        // Media progressiva sui primi n blocchi
        double progressive_average =
            progressive_sum
            / static_cast<double>(n);


        // Media progressiva dei quadrati
        double progressive_average2 =
            progressive_sum2
            / static_cast<double>(n);


        // Errore statistico progressivo
        double uncertainty =
            error(
                progressive_average,
                progressive_average2,
                n
            );


        // Scrittura dei risultati sul file
        energy_file
            << n << " "
            << block_average << " "
            << progressive_average << " "
            << uncertainty
            << '\n';
    }


    // --------------------------------------------------------
    // RISULTATI FINALI
    // --------------------------------------------------------

    double acceptance_rate =
        static_cast<double>(accepted)
        / static_cast<double>(attempted);


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


    // --------------------------------------------------------
    // STAMPA A TERMINALE
    // --------------------------------------------------------

    cout << endl;

    cout << "Variational Monte Carlo" << endl;

    cout << "-----------------------" << endl;

    cout << "mu      = "
         << mu
         << endl;

    cout << "sigma   = "
         << sigma
         << endl;

    cout << "delta   = "
         << delta
         << endl;

    cout << endl;


    cout << "<H> = "
         << final_energy
         << " +/- "
         << final_error
         << endl;


    cout << endl;


    cout << "Accettazione Metropolis = "
         << 100.0 * acceptance_rate
         << " %"
         << endl;


    // --------------------------------------------------------
    // CHIUSURA FILE E SALVATAGGIO DEL SEED
    // --------------------------------------------------------

    energy_file.close();

    sample_file.close();

    rnd.SaveSeed();


    return 0;
}