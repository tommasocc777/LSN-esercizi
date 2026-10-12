#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;


/*
    Exercise 02.1 - Monte Carlo integration

    Viene stimato numericamente l'integrale

        I = integral_0^1 [pi/2 * cos(pi*x/2)] dx = 1

    utilizzando due diversi metodi Monte Carlo:

    1. Campionamento uniforme nell'intervallo [0,1].
    2. Importance sampling con una distribuzione ottenuta
       dallo sviluppo di Taylor al primo ordine dell'integranda,
       centrato in x = 1/2.

    Le stime vengono analizzate tramite il metodo del data blocking.

    Il file integral.txt contiene, per ogni blocco:
    - numero di blocchi utilizzati;
    - stima progressiva con campionamento uniforme;
    - errore statistico della stima uniforme;
    - stima progressiva con importance sampling;
    - errore statistico della stima con importance sampling.
*/


// Funzione integranda
double integrand(double x) {

    const double pi = acos(-1.0);

    return (pi / 2.0) * cos(pi * x / 2.0);
}


// Distribuzione di probabilità usata per l'importance sampling
//
// Deriva dallo sviluppo di Taylor al primo ordine
// dell'integranda centrato in x = 1/2:
//
// p(x) = 1 - pi/2 * (x - 1/2)
double importance_pdf(double x) {

    const double pi = acos(-1.0);

    return 1.0 - (pi / 2.0) * (x - 0.5);
}


// Campionamento della distribuzione usata per l'importance sampling
double sample_importance(Random &rnd) {

    const double pi = acos(-1.0);

    // Viene estratto un numero uniforme in [0,1)
    double r = rnd.Rannyu();

    /*
        La funzione cumulativa è

            F(x) = (1 + pi/4)x - (pi/4)x^2.

        Viene utilizzata l'inversione della cumulativa
        per ottenere una variabile distribuita secondo p(x).
    */

    double a = 1.0 + pi / 4.0;

    return (2.0 / pi) *
           (a - sqrt(a * a - pi * r));
}


// Errore statistico della media progressiva
double error(double ave, double ave2, int n) {

    if (n == 0)
        return 0.0;

    return sqrt((ave2 - ave * ave) / n);
}


int main() {

    // ------------------------------------------------------------
    // Inizializzazione del generatore di numeri casuali
    // ------------------------------------------------------------

    Random rnd;

    int seed[4];
    int p1, p2;


    ifstream primes("Primes");

    if (primes.is_open()) {

        primes >> p1 >> p2;

    } else {

        cerr << "Errore: impossibile aprire il file Primes" << endl;
        return 1;
    }

    primes.close();


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

    } else {

        cerr << "Errore: impossibile aprire il file seed.in" << endl;
        return 1;
    }

    input.close();


    // ------------------------------------------------------------
    // Parametri della simulazione
    // ------------------------------------------------------------

    const int m = 1000000;        // Numero totale di estrazioni
    const int n_blocks = 100;     // Numero di blocchi
    const int l = m / n_blocks;   // Numero di estrazioni per blocco


    // ------------------------------------------------------------
    // Apertura del file di output
    // ------------------------------------------------------------

    ofstream output("integral.txt");

    if (!output.is_open()) {

        cerr << "Errore: impossibile creare integral.txt" << endl;
        return 1;
    }


    // ------------------------------------------------------------
    // Intestazione del file di output
    // ------------------------------------------------------------

    output << "# Exercise 02.1 - Monte Carlo integration" << endl;
    output << "# Stima di I = integral_0^1 [pi/2 * cos(pi*x/2)] dx = 1" << endl;
    output << "#" << endl;
    output << "# Metodo 1: campionamento uniforme in [0,1]" << endl;
    output << "# Metodo 2: importance sampling con distribuzione ottenuta" << endl;
    output << "#           dal Taylor al primo ordine centrato in x = 1/2" << endl;
    output << "#" << endl;
    output << "# M = " << m << " estrazioni totali" << endl;
    output << "# N = " << n_blocks << " blocchi" << endl;
    output << "# L = " << l << " estrazioni per blocco" << endl;
    output << "#" << endl;
    output << "# Colonne:" << endl;
    output << "# 1: numero di blocchi" << endl;
    output << "# 2: stima progressiva - campionamento uniforme" << endl;
    output << "# 3: errore statistico - campionamento uniforme" << endl;
    output << "# 4: stima progressiva - importance sampling" << endl;
    output << "# 5: errore statistico - importance sampling" << endl;
    output << "#" << endl;


    // ------------------------------------------------------------
    // Variabili per le medie progressive
    // ------------------------------------------------------------

    double sum_uniform = 0.0;
    double sum2_uniform = 0.0;

    double sum_importance = 0.0;
    double sum2_importance = 0.0;


    // ------------------------------------------------------------
    // Data blocking
    // ------------------------------------------------------------

    for (int i = 0; i < n_blocks; ++i) {

        double block_uniform = 0.0;
        double block_importance = 0.0;


        for (int j = 0; j < l; ++j) {

            // ----------------------------------------------------
            // Campionamento uniforme in [0,1]
            // ----------------------------------------------------

            double x_uniform = rnd.Rannyu();

            block_uniform += integrand(x_uniform);


            // ----------------------------------------------------
            // Importance sampling
            // ----------------------------------------------------

            double x_importance = sample_importance(rnd);

            /*
                L'integranda viene divisa per la densità
                utilizzata per il campionamento.
            */

            block_importance += integrand(x_importance)
                              / importance_pdf(x_importance);
        }


        // --------------------------------------------------------
        // Viene calcolata la media di ciascun blocco
        // --------------------------------------------------------

        block_uniform /= l;
        block_importance /= l;


        // --------------------------------------------------------
        // Vengono aggiornate le somme progressive
        // --------------------------------------------------------

        sum_uniform += block_uniform;
        sum2_uniform += block_uniform * block_uniform;

        sum_importance += block_importance;
        sum2_importance += block_importance * block_importance;


        // --------------------------------------------------------
        // Vengono calcolate le medie progressive
        // --------------------------------------------------------

        double ave_uniform = sum_uniform / (i + 1);
        double ave2_uniform = sum2_uniform / (i + 1);

        double ave_importance = sum_importance / (i + 1);
        double ave2_importance = sum2_importance / (i + 1);


        // --------------------------------------------------------
        // Vengono calcolati gli errori statistici progressivi
        // --------------------------------------------------------

        double err_uniform =
            error(ave_uniform, ave2_uniform, i);

        double err_importance =
            error(ave_importance, ave2_importance, i);


        // --------------------------------------------------------
        // Salvataggio dei risultati
        // --------------------------------------------------------

        output << i + 1 << " "
               << ave_uniform << " "
               << err_uniform << " "
               << ave_importance << " "
               << err_importance << endl;
    }


    output.close();

    rnd.SaveSeed();

    return 0;
}