/***********************************************************************
 Esercizio 01.2 - Teorema del Limite Centrale

 Il programma genera 10000 realizzazioni della variabile

       S_N = (1/N) * somma_i x_i

 per N = 1, 2, 10, 100 e per tre diverse distribuzioni:

 1) distribuzione uniforme in [0,1);
 2) distribuzione esponenziale con lambda = 1;
 3) distribuzione di Cauchy-Lorentz con mu = 0 e Gamma = 1.

 File di output:
    Standard_dice.txt
    Exponential_dice.txt
    Lorentzian_dice.txt

 Ogni file contiene quattro colonne, corrispondenti a:
    N = 1, 2, 10, 100
************************************************************************/

#include <iostream>     
#include <fstream>      
#include <string>       
#include "random.h"     
using namespace std;


int main() {

   // =================================================================
   // INIZIALIZZAZIONE DEL GENERATORE PSEUDO-CASUALE
   // =================================================================

   Random rnd;                       // Generatore pseudo-casuale
   int seed[4];                      // Quattro interi che definiscono il seed
   int primo_1, primo_2;             // Numeri primi utilizzati dal generatore

   ifstream file_primes("Primes");   // Apro il file contenente i numeri primi

   if (file_primes.is_open()) {

      file_primes >> primo_1 >> primo_2;   // Leggo i due numeri primi
   }
   else {

      cerr << "ERRORE: impossibile aprire il file Primes" << endl;
      return 1;
   }

   file_primes.close();   // Chiudo il file Primes


   ifstream file_seed("seed.in");   // Apro il file contenente il seed

   string proprieta;                // Stringa utilizzata per cercare RANDOMSEED
   bool seed_trovato = false;       // Controlla che il seed sia stato trovato

   if (file_seed.is_open()) {

      while (file_seed >> proprieta) {   // Leggo il file parola per parola

         if (proprieta == "RANDOMSEED") {

            file_seed >> seed[0]
                      >> seed[1]
                      >> seed[2]
                      >> seed[3];         // Leggo i quattro valori del seed

            rnd.SetRandom(seed, primo_1, primo_2);   // Inizializzo il generatore

            seed_trovato = true;   // Segnalo che il seed e' stato trovato
         }
      }

      file_seed.close();   // Chiudo seed.in
   }
   else {

      cerr << "ERRORE: impossibile aprire il file seed.in" << endl;
      return 1;
   }

   if (!seed_trovato) {

      cerr << "ERRORE: RANDOMSEED non trovato in seed.in" << endl;
      return 1;
   }


   // =================================================================
   // PARAMETRI DELLA SIMULAZIONE
   // =================================================================

   const int numero_realizzazioni = 10000;   // Numero di valori di S_N da generare

   const int numero_n = 4;                   // Numero di diversi valori di N

   const int valori_n[numero_n] =
      {1, 2, 10, 100};                       // Valori di N richiesti dalla consegna

   const double lambda = 1.0;               // Parametro della distribuzione esponenziale

   const double mu = 0.0;                   // Centro della distribuzione Lorentziana

   const double gamma = 1.0;                // Larghezza della distribuzione Lorentziana


   // =================================================================
   // APERTURA DEI FILE DI OUTPUT
   // =================================================================

   ofstream file_standard("Standard_dice.txt");       // Distribuzione uniforme
   ofstream file_exponential("Exponential_dice.txt"); // Distribuzione esponenziale
   ofstream file_lorentzian("Lorentzian_dice.txt");   // Distribuzione Lorentziana

   if (!file_standard.is_open() ||
       !file_exponential.is_open() ||
       !file_lorentzian.is_open()) {

      cerr << "ERRORE: impossibile creare uno o piu' file di output" << endl;
      return 1;
   }


   // Header dei file: le righe con # verranno ignorate da numpy.loadtxt()

   file_standard << "# Esercizio 01.2 - Distribuzione uniforme" << endl;
   file_standard << "# Colonna 1: N = 1" << endl;
   file_standard << "# Colonna 2: N = 2" << endl;
   file_standard << "# Colonna 3: N = 10" << endl;
   file_standard << "# Colonna 4: N = 100" << endl;

   file_exponential << "# Esercizio 01.2 - Distribuzione esponenziale" << endl;
   file_exponential << "# lambda = " << lambda << endl;
   file_exponential << "# Colonna 1: N = 1" << endl;
   file_exponential << "# Colonna 2: N = 2" << endl;
   file_exponential << "# Colonna 3: N = 10" << endl;
   file_exponential << "# Colonna 4: N = 100" << endl;

   file_lorentzian << "# Esercizio 01.2 - Distribuzione Cauchy-Lorentz" << endl;
   file_lorentzian << "# mu = " << mu
                   << ", gamma = " << gamma << endl;
   file_lorentzian << "# Colonna 1: N = 1" << endl;
   file_lorentzian << "# Colonna 2: N = 2" << endl;
   file_lorentzian << "# Colonna 3: N = 10" << endl;
   file_lorentzian << "# Colonna 4: N = 100" << endl;


   // =================================================================
   // GENERAZIONE DELLE 10000 REALIZZAZIONI DI S_N
   // =================================================================

   for (int realizzazione = 0;
        realizzazione < numero_realizzazioni;
        realizzazione++) {

      // Per ogni riga calcoliamo S_N per N = 1, 2, 10, 100
      for (int k = 0; k < numero_n; k++) {

         int n = valori_n[k];   // Numero di variabili da mediare

         double somma_standard = 0.0;       // Somma dei numeri uniformi
         double somma_exponential = 0.0;    // Somma dei numeri esponenziali
         double somma_lorentzian = 0.0;     // Somma dei numeri Lorentziani


         // Generiamo N numeri per ciascuna distribuzione
         for (int i = 0; i < n; i++) {

            somma_standard += rnd.Rannyu();
            // Variabile uniforme in [0,1)

            somma_exponential += rnd.Exponential(lambda);
            // Variabile esponenziale con parametro lambda

            somma_lorentzian += rnd.Lorentz(mu, gamma);
            // Variabile di Cauchy-Lorentz con parametri mu e gamma
         }


         double s_n_standard =
            somma_standard / n;          // S_N per la distribuzione uniforme

         double s_n_exponential =
            somma_exponential / n;       // S_N per la distribuzione esponenziale

         double s_n_lorentzian =
            somma_lorentzian / n;        // S_N per la distribuzione Lorentziana


         file_standard << s_n_standard;           // Scrivo S_N uniforme
         file_exponential << s_n_exponential;     // Scrivo S_N esponenziale
         file_lorentzian << s_n_lorentzian;       // Scrivo S_N Lorentziano


         // Se non siamo all'ultima colonna aggiungo uno spazio
         if (k < numero_n - 1) {

            file_standard << " ";
            file_exponential << " ";
            file_lorentzian << " ";
         }
      }


      // Fine della realizzazione: vado a capo nei tre file
      file_standard << endl;
      file_exponential << endl;
      file_lorentzian << endl;
   }


   // =================================================================
   // CHIUSURA DEI FILE E DEL PROGRAMMA
   // =================================================================

   file_standard.close();       // Chiudo il file della distribuzione uniforme
   file_exponential.close();    // Chiudo il file della distribuzione esponenziale
   file_lorentzian.close();     // Chiudo il file della distribuzione Lorentziana

   rnd.SaveSeed();              // Salvo lo stato finale del generatore

   return 0;                    // Fine corretta del programma
}