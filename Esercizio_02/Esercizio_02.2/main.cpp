#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"

using namespace std;


/*
    Exercise 02.2 - Random Walks

    PARTE PRINCIPALE

    Vengono simulati due random walk tridimensionali:

    1. Random walk su un reticolo cubico con passo a = 1.
    2. Random walk nel continuo con passo a = 1.

    Per ogni numero di step N viene calcolata

        sqrt(<|r_N|^2>)

    utilizzando il metodo del data blocking.

    Il file random_walk.txt contiene:
    - numero di step;
    - RMS del random walk sul reticolo;
    - errore statistico sul reticolo;
    - RMS del random walk nel continuo;
    - errore statistico nel continuo.


    PARTE EXTRA

    Viene studiata la probabilita' cumulativa di ritorno
    all'origine per random walk su reticolo in:

    - 1 dimensione;
    - 2 dimensioni;
    - 3 dimensioni.

    Per ogni numero di step N viene calcolata la probabilita'
    che il walker sia ritornato almeno una volta all'origine
    entro lo step N.

    Il file return_probability.txt contiene:
    - numero di step;
    - probabilita' di ritorno in 1D e relativo errore;
    - probabilita' di ritorno in 2D e relativo errore;
    - probabilita' di ritorno in 3D e relativo errore.
*/


// ------------------------------------------------------------
// Errore statistico della media
// ------------------------------------------------------------

double error(double ave, double ave2, int n) {

    if (n == 0)
        return 0.0;

    double variance = ave2 - ave * ave;

    // Viene evitato un eventuale piccolo valore negativo
    // dovuto alla precisione numerica.
    if (variance < 0.0)
        variance = 0.0;

    return sqrt(variance / n);
}


// ------------------------------------------------------------
// Passo su un reticolo in d dimensioni
// ------------------------------------------------------------

void lattice_step(Random &rnd,
                  int dimension,
                  int &x,
                  int &y,
                  int &z) {

    /*
        In d dimensioni sono disponibili 2*d direzioni.

        d = 1:
            +x, -x

        d = 2:
            +x, -x, +y, -y

        d = 3:
            +x, -x, +y, -y, +z, -z
    */

    int direction =
        static_cast<int>(2.0 * dimension * rnd.Rannyu());


    if (direction == 0)
        x += 1;

    else if (direction == 1)
        x -= 1;

    else if (direction == 2)
        y += 1;

    else if (direction == 3)
        y -= 1;

    else if (direction == 4)
        z += 1;

    else
        z -= 1;
}


// ------------------------------------------------------------
// Passo nel continuo tridimensionale
// ------------------------------------------------------------

void continuum_step(Random &rnd,
                    double &x,
                    double &y,
                    double &z,
                    double a) {

    const double pi = acos(-1.0);

    /*
        Per campionare uniformemente una direzione nello spazio:

            phi = 2*pi*r1

        mentre theta deve essere distribuito secondo

            p(theta) = 1/2 * sin(theta).

        Invertendo la cumulativa:

            theta = arccos(1 - 2*r2),

        con r1 e r2 uniformi in [0,1).
    */

    double phi =
        2.0 * pi * rnd.Rannyu();

    double r =
        rnd.Rannyu();

    double theta =
        acos(1.0 - 2.0 * r);


    // Componenti cartesiane del passo
    x += a * sin(theta) * cos(phi);
    y += a * sin(theta) * sin(phi);
    z += a * cos(theta);
}


// ============================================================
// PARTE EXTRA
// Probabilita' cumulativa di ritorno all'origine
// ============================================================

void return_probability(Random &rnd,
                        int m,
                        int n_blocks,
                        int n_steps,
                        const string &filename) {

    const int l = m / n_blocks;


    /*
        Le tre righe corrispondono rispettivamente a:

            index 0 -> 1D
            index 1 -> 2D
            index 2 -> 3D
    */

    vector<vector<double>> sum_return(
        3,
        vector<double>(n_steps + 1, 0.0)
    );

    vector<vector<double>> sum2_return(
        3,
        vector<double>(n_steps + 1, 0.0)
    );


    // --------------------------------------------------------
    // Ciclo sui blocchi
    // --------------------------------------------------------

    for (int iblock = 0; iblock < n_blocks; ++iblock) {


        /*
            new_returns[d][N] contiene il numero di walker
            che ritornano all'origine per la prima volta
            allo step N.

            Questa informazione viene utilizzata solamente
            internamente per costruire la probabilita'
            cumulativa di ritorno.

            La distribuzione del primo ritorno non viene
            salvata né analizzata.
        */

        vector<vector<int>> new_returns(
            3,
            vector<int>(n_steps + 1, 0)
        );


        // ----------------------------------------------------
        // Simulazioni in 1D, 2D e 3D
        // ----------------------------------------------------

        for (int idim = 0; idim < 3; ++idim) {

            int dimension = idim + 1;


            // -----------------------------------------------
            // Random walk del blocco corrente
            // -----------------------------------------------

            for (int iwalk = 0; iwalk < l; ++iwalk) {

                int x = 0;
                int y = 0;
                int z = 0;


                /*
                    Lo stato iniziale N = 0 non viene
                    considerato come un ritorno.
                */

                for (int istep = 1;
                     istep <= n_steps;
                     ++istep) {

                    lattice_step(
                        rnd,
                        dimension,
                        x,
                        y,
                        z
                    );


                    // ---------------------------------------
                    // Controllo del ritorno all'origine
                    // ---------------------------------------

                    if (x == 0 &&
                        y == 0 &&
                        z == 0) {

                        new_returns[idim][istep] += 1;

                        /*
                            Una volta avvenuto un ritorno,
                            il walker appartiene alla
                            probabilita' cumulativa per tutti
                            gli step successivi.

                            Non e' quindi necessario
                            continuare il random walk.
                        */

                        break;
                    }
                }
            }


            // -----------------------------------------------
            // Probabilita' cumulativa nel blocco corrente
            // -----------------------------------------------

            int cumulative_returns = 0;


            for (int istep = 1;
                 istep <= n_steps;
                 ++istep) {

                cumulative_returns +=
                    new_returns[idim][istep];


                double p_return =
                    static_cast<double>(
                        cumulative_returns
                    ) / l;


                // Accumulo per il data blocking
                sum_return[idim][istep] +=
                    p_return;

                sum2_return[idim][istep] +=
                    p_return * p_return;
            }
        }
    }


    // --------------------------------------------------------
    // Apertura del file di output
    // --------------------------------------------------------

    ofstream output(filename);

    if (!output.is_open()) {

        cerr << "Errore: impossibile creare "
             << filename << endl;

        return;
    }


    // --------------------------------------------------------
    // Intestazione del file
    // --------------------------------------------------------

    output << "# Exercise 02.2 - EXTRA" << endl;
    output << "# Cumulative return probability" << endl;
    output << "#" << endl;

    output << "# Viene studiata la probabilita' che un random walk"
           << endl;

    output << "# sia ritornato almeno una volta all'origine"
           << endl;

    output << "# entro lo step N." << endl;
    output << "#" << endl;

    output << "# Vengono confrontati random walk su reticolo"
           << endl;

    output << "# in 1D, 2D e 3D." << endl;
    output << "#" << endl;

    output << "# Lo stato iniziale N = 0 non viene considerato"
           << endl;

    output << "# come ritorno all'origine." << endl;
    output << "#" << endl;

    output << "# M = "
           << m
           << " random walk per dimensione"
           << endl;

    output << "# N_blocks = "
           << n_blocks
           << endl;

    output << "# L = "
           << l
           << " random walk per blocco"
           << endl;

    output << "# Step massimi = "
           << n_steps
           << endl;

    output << "#" << endl;

    output << "# Colonne:" << endl;

    output << "# 1: numero di step" << endl;

    output << "# 2: probabilita' cumulativa di ritorno - 1D"
           << endl;

    output << "# 3: errore statistico - 1D"
           << endl;

    output << "# 4: probabilita' cumulativa di ritorno - 2D"
           << endl;

    output << "# 5: errore statistico - 2D"
           << endl;

    output << "# 6: probabilita' cumulativa di ritorno - 3D"
           << endl;

    output << "# 7: errore statistico - 3D"
           << endl;

    output << "#" << endl;


    // --------------------------------------------------------
    // N = 0
    // --------------------------------------------------------

    output << 0 << " "
           << 0.0 << " "
           << 0.0 << " "
           << 0.0 << " "
           << 0.0 << " "
           << 0.0 << " "
           << 0.0 << endl;


    // --------------------------------------------------------
    // Media finale ed errore sui blocchi
    // --------------------------------------------------------

    for (int istep = 1;
         istep <= n_steps;
         ++istep) {

        double results[3];
        double errors[3];


        for (int idim = 0;
             idim < 3;
             ++idim) {

            double ave =
                sum_return[idim][istep]
                / n_blocks;

            double ave2 =
                sum2_return[idim][istep]
                / n_blocks;


            results[idim] = ave;

            errors[idim] =
                error(
                    ave,
                    ave2,
                    n_blocks - 1
                );
        }


        output << istep << " "
               << results[0] << " "
               << errors[0] << " "
               << results[1] << " "
               << errors[1] << " "
               << results[2] << " "
               << errors[2] << endl;
    }


    output.close();
}


// ============================================================
// MAIN
// ============================================================

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

        cerr << "Errore: impossibile aprire il file Primes"
             << endl;

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

                rnd.SetRandom(
                    seed,
                    p1,
                    p2
                );
            }
        }

    } else {

        cerr << "Errore: impossibile aprire il file seed.in"
             << endl;

        return 1;
    }

    input.close();


    // ============================================================
    // PARTE PRINCIPALE
    // ============================================================

    const int m = 10000;
    const int n_blocks = 100;

    const int l =
        m / n_blocks;

    const int n_steps = 100;

    const double a = 1.0;


    // ------------------------------------------------------------
    // Apertura del file
    // ------------------------------------------------------------

    ofstream output(
        "random_walk.txt"
    );


    if (!output.is_open()) {

        cerr << "Errore: impossibile creare random_walk.txt"
             << endl;

        return 1;
    }


    // ------------------------------------------------------------
    // Intestazione
    // ------------------------------------------------------------

    output << "# Exercise 02.2 - 3D Random Walks"
           << endl;

    output << "#" << endl;

    output << "# Vengono simulati due random walk tridimensionali:"
           << endl;

    output << "# 1. reticolo cubico con passo a = 1"
           << endl;

    output << "# 2. continuo con direzione campionata uniformemente"
           << endl;

    output << "#    sull'angolo solido"
           << endl;

    output << "#" << endl;

    output << "# Nel continuo:"
           << endl;

    output << "# phi   = 2*pi*r1"
           << endl;

    output << "# theta = arccos(1 - 2*r2)"
           << endl;

    output << "#" << endl;

    output << "# Per ogni step viene calcolata"
           << endl;

    output << "# sqrt(<|r|^2>) tramite data blocking."
           << endl;

    output << "#" << endl;

    output << "# M = "
           << m
           << endl;

    output << "# N_blocks = "
           << n_blocks
           << endl;

    output << "# L = "
           << l
           << endl;

    output << "# Step massimi = "
           << n_steps
           << endl;

    output << "# a = "
           << a
           << endl;

    output << "#" << endl;

    output << "# Colonne:"
           << endl;

    output << "# 1: numero di step"
           << endl;

    output << "# 2: sqrt(<r^2>) - reticolo cubico"
           << endl;

    output << "# 3: errore statistico - reticolo cubico"
           << endl;

    output << "# 4: sqrt(<r^2>) - continuo"
           << endl;

    output << "# 5: errore statistico - continuo"
           << endl;

    output << "#" << endl;


    // ------------------------------------------------------------
    // Variabili per il data blocking
    // ------------------------------------------------------------

    vector<double> sum_lattice(
        n_steps + 1,
        0.0
    );

    vector<double> sum2_lattice(
        n_steps + 1,
        0.0
    );

    vector<double> sum_continuum(
        n_steps + 1,
        0.0
    );

    vector<double> sum2_continuum(
        n_steps + 1,
        0.0
    );


    // ------------------------------------------------------------
    // Ciclo sui blocchi
    // ------------------------------------------------------------

    for (int iblock = 0;
         iblock < n_blocks;
         ++iblock) {

        vector<double> block_r2_lattice(
            n_steps + 1,
            0.0
        );

        vector<double> block_r2_continuum(
            n_steps + 1,
            0.0
        );


        for (int iwalk = 0;
             iwalk < l;
             ++iwalk) {

            // Reticolo
            int x_lattice = 0;
            int y_lattice = 0;
            int z_lattice = 0;


            // Continuo
            double x_continuum = 0.0;
            double y_continuum = 0.0;
            double z_continuum = 0.0;


            for (int istep = 1;
                 istep <= n_steps;
                 ++istep) {


                // -------------------------------------------
                // Random walk sul reticolo 3D
                // -------------------------------------------

                lattice_step(
                    rnd,
                    3,
                    x_lattice,
                    y_lattice,
                    z_lattice
                );


                // -------------------------------------------
                // Random walk nel continuo 3D
                // -------------------------------------------

                continuum_step(
                    rnd,
                    x_continuum,
                    y_continuum,
                    z_continuum,
                    a
                );


                // -------------------------------------------
                // Distanze quadratiche
                // -------------------------------------------

                double r2_lattice =
                    a * a *
                    (
                        x_lattice * x_lattice +
                        y_lattice * y_lattice +
                        z_lattice * z_lattice
                    );


                double r2_continuum =
                    x_continuum * x_continuum +
                    y_continuum * y_continuum +
                    z_continuum * z_continuum;


                block_r2_lattice[istep] +=
                    r2_lattice;

                block_r2_continuum[istep] +=
                    r2_continuum;
            }
        }


        // --------------------------------------------------------
        // Media nel blocco
        // --------------------------------------------------------

        for (int istep = 0;
             istep <= n_steps;
             ++istep) {

            double mean_r2_lattice =
                block_r2_lattice[istep]
                / l;

            double mean_r2_continuum =
                block_r2_continuum[istep]
                / l;


            double rms_lattice =
                sqrt(mean_r2_lattice);

            double rms_continuum =
                sqrt(mean_r2_continuum);


            sum_lattice[istep] +=
                rms_lattice;

            sum2_lattice[istep] +=
                rms_lattice * rms_lattice;


            sum_continuum[istep] +=
                rms_continuum;

            sum2_continuum[istep] +=
                rms_continuum * rms_continuum;
        }
    }


    // ------------------------------------------------------------
    // Media finale ed errore
    // ------------------------------------------------------------

    for (int istep = 0;
         istep <= n_steps;
         ++istep) {


        double ave_lattice =
            sum_lattice[istep]
            / n_blocks;

        double ave2_lattice =
            sum2_lattice[istep]
            / n_blocks;


        double ave_continuum =
            sum_continuum[istep]
            / n_blocks;

        double ave2_continuum =
            sum2_continuum[istep]
            / n_blocks;


        double err_lattice =
            error(
                ave_lattice,
                ave2_lattice,
                n_blocks - 1
            );


        double err_continuum =
            error(
                ave_continuum,
                ave2_continuum,
                n_blocks - 1
            );


        output << istep << " "
               << ave_lattice << " "
               << err_lattice << " "
               << ave_continuum << " "
               << err_continuum
               << endl;
    }


    output.close();


    // ============================================================
    // PARTE EXTRA
    // Probabilita' cumulativa di ritorno
    // ============================================================

    const int m_extra = 50000;

    const int n_blocks_extra = 100;

    const int n_steps_extra = 2000;


    return_probability(
        rnd,
        m_extra,
        n_blocks_extra,
        n_steps_extra,
        "return_probability.txt"
    );


    // ------------------------------------------------------------
    // Salvataggio del seed finale
    // ------------------------------------------------------------

    rnd.SaveSeed();

    return 0;
}