/*
===============================================================================
ESERCIZIO 03.1 - PRICING DI OPZIONI EUROPEE CON MONTE CARLO
===============================================================================

Obiettivo:
stimare il prezzo al tempo t = 0 di una European Call e di una European Put
assumendo che il prezzo dell'asset segua un Geometric Brownian Motion:

    dS(t) = r S(t) dt + sigma S(t) dW(t)

I payoff a scadenza sono:

    Call: max[0, S(T) - K]
    Put : max[0, K - S(T)]

e vengono scontati al tempo t = 0 tramite il fattore exp(-rT).

La simulazione viene eseguita in due modi:

1) campionamento diretto di S(T);
2) campionamento della traiettoria in 100 intervalli temporali.

Le stime vengono analizzate tramite data blocking e confrontate con
la soluzione analitica di Black-Scholes.

Ogni file di output contiene:

    colonna 1 -> media del singolo blocco
    colonna 2 -> media progressiva
    colonna 3 -> errore progressivo

===============================================================================
*/


#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
#include <filesystem>
#include <algorithm>

#include "random.h"

using namespace std;


// Funzione di distribuzione cumulativa della normale standard
double normal_cdf(double x) {

    return 0.5 * (1.0 + erf(x / sqrt(2.0)));
}


// ============================================================================
// SCRITTURA DEI RISULTATI DEL DATA BLOCKING
// ============================================================================

void write_blocking(const string& filename,
                    const vector<double>& block_values,
                    const string& tipo_opzione,
                    const string& metodo_campionamento) {

    ofstream output(filename);

    if (!output.is_open()) {

        cerr << "Errore: impossibile aprire il file "
             << filename << endl;

        return;
    }


    // ------------------------------------------------------------------------
    // Intestazione del file
    // ------------------------------------------------------------------------

    output << "# ================================================================\n";
    output << "# Esercizio 03.1 - Pricing di opzioni europee con Monte Carlo\n";
    output << "#\n";

    output << "# Tipo di opzione: " << tipo_opzione << "\n";

    output << "# Metodo di campionamento: "
           << metodo_campionamento << "\n";

    output << "#\n";

    output << "# La quantita' riportata rappresenta il payoff scontato\n";
    output << "# dell'opzione, stimato mediante il metodo del data blocking.\n";

    output << "#\n";

    output << "# Colonne:\n";
    output << "# 1 -> media_blocco       : media del singolo blocco\n";
    output << "# 2 -> media_progressiva  : media calcolata usando i blocchi disponibili\n";
    output << "# 3 -> errore_progressivo : incertezza statistica sulla media progressiva\n";

    output << "# ================================================================\n";

    output << "# media_blocco media_progressiva errore_progressivo\n";


    // ------------------------------------------------------------------------
    // Variabili per il calcolo progressivo
    // ------------------------------------------------------------------------

    double sum = 0.0;          // Somma delle medie dei blocchi
    double sum2 = 0.0;         // Somma dei quadrati delle medie dei blocchi


    for (size_t i = 0; i < block_values.size(); ++i) {

        // Media del blocco corrente
        double value = block_values[i];


        // Aggiornamento delle somme progressive
        sum += value;

        sum2 += value * value;


        // Numero di blocchi utilizzati fino a questo punto
        double n =
            static_cast<double>(i + 1);


        // Media progressiva
        double progressive_mean =
            sum / n;


        // Media progressiva dei quadrati
        double progressive_mean2 =
            sum2 / n;


        // Errore statistico progressivo
        double error = 0.0;


        // Con un solo blocco non possiamo stimare l'errore
        if (i > 0) {

            double variance =
                progressive_mean2
                -
                progressive_mean * progressive_mean;


            // Protezione da piccoli valori negativi
            // dovuti alla precisione numerica
            if (variance < 0.0) {

                variance = 0.0;
            }


            error =
                sqrt(
                    variance
                    /
                    static_cast<double>(i)
                );
        }


        // Scrittura dei risultati
        output << scientific << setprecision(10)
               << value << " "
               << progressive_mean << " "
               << error << "\n";
    }


    output.close();
}


// ============================================================================
// INIZIALIZZAZIONE DEL GENERATORE DI NUMERI PSEUDO-CASUALI
// ============================================================================

void initialize_random(Random& rnd) {

    int seed[4];

    int p1;
    int p2;


    // Lettura dei numeri primi
    ifstream primes("Primes");


    if (primes.is_open()) {

        primes >> p1 >> p2;

        primes.close();

    } else {

        cerr << "Errore: impossibile aprire il file Primes" << endl;

        exit(1);
    }


    // Lettura del seed
    ifstream input("seed.in");

    string property;


    if (input.is_open()) {

        while (input >> property) {

            if (property == "RANDOMSEED") {

                input >> seed[0]
                      >> seed[1]
                      >> seed[2]
                      >> seed[3];


                rnd.SetRandom(seed, p1, p2);
            }
        }


        input.close();

    } else {

        cerr << "Errore: impossibile aprire il file seed.in" << endl;

        exit(1);
    }
}


int main() {

    // =========================================================================
    // PARAMETRI DEL PROBLEMA
    // =========================================================================

    const double S0 = 100.0;       // Prezzo iniziale dell'asset

    const double T = 1.0;          // Tempo di scadenza

    const double K = 100.0;        // Strike price

    const double r = 0.1;          // Tasso risk-free

    const double sigma = 0.25;     // Volatilita'


    // =========================================================================
    // PARAMETRI MONTE CARLO E DATA BLOCKING
    // =========================================================================

    const int M = 1000000;          // Numero totale di simulazioni

    const int N_blocks = 1000;      // Numero di blocchi

    const int L =
        M / N_blocks;              // Simulazioni per blocco


    // =========================================================================
    // PARAMETRI DEL CAMPIONAMENTO DISCRETIZZATO
    // =========================================================================

    const int N_steps = 100;       // Numero di intervalli temporali

    const double dt =
        T / static_cast<double>(N_steps);


    // =========================================================================
    // FATTORE DI SCONTO
    // =========================================================================

    // Riporta il payoff dal tempo T al tempo t = 0
    const double discount =
        exp(-r * T);


    // Controllo sulla divisione in blocchi
    if (M % N_blocks != 0) {

        cerr << "Errore: M deve essere divisibile per N_blocks" << endl;

        return 1;
    }


    // =========================================================================
    // GENERATORE DI NUMERI PSEUDO-CASUALI
    // =========================================================================

    Random rnd;

    initialize_random(rnd);


    // =========================================================================
    // VETTORI DELLE MEDIE DI BLOCCO
    // =========================================================================

    vector<double> call_direct(N_blocks);

    vector<double> put_direct(N_blocks);

    vector<double> call_discrete(N_blocks);

    vector<double> put_discrete(N_blocks);


    // =========================================================================
    // SIMULAZIONE MONTE CARLO
    // =========================================================================

    for (int iblock = 0; iblock < N_blocks; ++iblock) {

        // Accumulatori del blocco corrente
        double sum_call_direct = 0.0;

        double sum_put_direct = 0.0;

        double sum_call_discrete = 0.0;

        double sum_put_discrete = 0.0;


        // Simulazioni appartenenti al blocco corrente
        for (int i = 0; i < L; ++i) {


            // =================================================================
            // 1) CAMPIONAMENTO DIRETTO DI S(T)
            // =================================================================

            // Variabile normale standard Z ~ N(0,1)
            double Z =
                rnd.Gauss(0.0, 1.0);


            // Soluzione esatta del GBM direttamente al tempo T
            double S_T_direct =
                S0
                *
                exp(
                    (r - 0.5 * sigma * sigma) * T
                    +
                    sigma * Z * sqrt(T)
                );


            // Payoff della Call
            double payoff_call_direct =
                max(
                    0.0,
                    S_T_direct - K
                );


            // Payoff della Put
            double payoff_put_direct =
                max(
                    0.0,
                    K - S_T_direct
                );


            // Sconto dei payoff al tempo t = 0
            sum_call_direct +=
                discount * payoff_call_direct;


            sum_put_direct +=
                discount * payoff_put_direct;



            // =================================================================
            // 2) CAMPIONAMENTO DELLA TRAIETTORIA DEL GBM
            // =================================================================

            // Partiamo dal prezzo iniziale
            double S = S0;


            // Evoluzione:
            //
            // S(0) -> S(t1) -> ... -> S(T)

            for (int istep = 0;
                 istep < N_steps;
                 ++istep) {


                // Nuova variabile normale indipendente
                double Z_step =
                    rnd.Gauss(0.0, 1.0);


                // Evoluzione esatta del GBM sul singolo intervallo dt
                S *=
                    exp(
                        (r - 0.5 * sigma * sigma) * dt
                        +
                        sigma * Z_step * sqrt(dt)
                    );
            }


            // Dopo N_steps abbiamo raggiunto il tempo T
            double S_T_discrete = S;


            // Payoff della Call
            double payoff_call_discrete =
                max(
                    0.0,
                    S_T_discrete - K
                );


            // Payoff della Put
            double payoff_put_discrete =
                max(
                    0.0,
                    K - S_T_discrete
                );


            // Sconto dei payoff al tempo t = 0
            sum_call_discrete +=
                discount * payoff_call_discrete;


            sum_put_discrete +=
                discount * payoff_put_discrete;
        }


        // =====================================================================
        // MEDIE DEL BLOCCO
        // =====================================================================

        call_direct[iblock] =
            sum_call_direct
            /
            static_cast<double>(L);


        put_direct[iblock] =
            sum_put_direct
            /
            static_cast<double>(L);


        call_discrete[iblock] =
            sum_call_discrete
            /
            static_cast<double>(L);


        put_discrete[iblock] =
            sum_put_discrete
            /
            static_cast<double>(L);
    }


    // =========================================================================
    // CREAZIONE DELLA CARTELLA OUTPUT
    // =========================================================================

    filesystem::create_directories("output");


    // =========================================================================
    // SCRITTURA DEI RISULTATI
    // =========================================================================


    // Call - campionamento diretto
    write_blocking(
        "output/call_direct.txt",
        call_direct,
        "Call europea",
        "Campionamento diretto del prezzo finale S(T)"
    );


    // Put - campionamento diretto
    write_blocking(
        "output/put_direct.txt",
        put_direct,
        "Put europea",
        "Campionamento diretto del prezzo finale S(T)"
    );


    // Call - campionamento della traiettoria
    write_blocking(
        "output/call_discrete.txt",
        call_discrete,
        "Call europea",
        "Campionamento della traiettoria del GBM in 100 intervalli temporali"
    );


    // Put - campionamento della traiettoria
    write_blocking(
        "output/put_discrete.txt",
        put_discrete,
        "Put europea",
        "Campionamento della traiettoria del GBM in 100 intervalli temporali"
    );


    // =========================================================================
    // SOLUZIONE ANALITICA DI BLACK-SCHOLES
    // =========================================================================

    double d1 =
        (
            log(S0 / K)
            +
            (r + 0.5 * sigma * sigma) * T
        )
        /
        (
            sigma * sqrt(T)
        );


    double d2 =
        d1
        -
        sigma * sqrt(T);


    // Prezzo analitico della Call
    double call_BS =
        S0 * normal_cdf(d1)
        -
        K * exp(-r * T) * normal_cdf(d2);


    // Prezzo analitico della Put
    double put_BS =
        K * exp(-r * T) * normal_cdf(-d2)
        -
        S0 * normal_cdf(-d1);


    // =========================================================================
    // STAMPA A TERMINALE
    // =========================================================================

    cout << fixed << setprecision(6);


    cout << "\n";
    cout << "=========================================\n";
    cout << "  ESERCIZIO 03.1 - PRICING OPZIONI\n";
    cout << "=========================================\n\n";


    cout << "PARAMETRI:\n\n";

    cout << "S(0)       = " << S0 << "\n";

    cout << "K          = " << K << "\n";

    cout << "T          = " << T << "\n";

    cout << "r          = " << r << "\n";

    cout << "sigma      = " << sigma << "\n\n";


    cout << "SIMULAZIONE:\n\n";

    cout << "M          = " << M << "\n";

    cout << "N_blocks   = " << N_blocks << "\n";

    cout << "L          = " << L << "\n";

    cout << "N_steps    = " << N_steps << "\n\n";


    cout << "SOLUZIONE ANALITICA DI BLACK-SCHOLES:\n\n";

    cout << "Call = " << call_BS << "\n";

    cout << "Put  = " << put_BS << "\n\n";


    cout << "Risultati Monte Carlo salvati nella cartella output/\n";


    // =========================================================================
    // SALVATAGGIO DEL SEED
    // =========================================================================

    rnd.SaveSeed();


    return 0;
}