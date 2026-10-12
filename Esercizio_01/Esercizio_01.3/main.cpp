/***********************************************************************
 Esercizio 01.3 - Esperimento dell'ago di Buffon

 Simuliamo il lancio casuale di un ago di lunghezza L su un piano
 formato da linee parallele poste a distanza d = 2L.

 La probabilita' che l'ago intersechi una linea e':

       P = 2L / (pi d)

 da cui possiamo stimare pi come:

       pi = 2L N_thr / (d N_hit)

 L'orientazione dell'ago viene generata senza utilizzare pi:
 scegliamo un punto casuale all'interno del disco unitario e
 normalizziamo il vettore ottenuto.

 Vengono eseguiti 10^8 lanci totali, suddivisi in blocchi
 da 10^4 lanci ciascuno.

 File di output:
    Pi_Buffon.txt
************************************************************************/

#include <iostream>     
#include <fstream>     
#include <string>       
#include <vector>       
#include <cmath>        
#include <algorithm>    
#include "random.h"     

using namespace std;


// Calcola l'errore statistico della media progressiva
double errore_progressivo(double media,
                          double media_quadrata,
                          int n_blocchi) {

   if (n_blocchi <= 1) return 0.0;   // Con un solo blocco l'errore non e' stimabile

   double varianza =
      (media_quadrata - media * media) / (n_blocchi - 1);   // Varianza della media

   if (varianza < 0.0) varianza = 0.0;   // Protezione da piccoli errori numerici

   return sqrt(varianza);   // Incertezza statistica
}


// Genera la componente verticale di una direzione casuale nel piano
double direzione_y_casuale(Random& rnd) {

   double x_dir;     // Componente x del vettore casuale
   double y_dir;     // Componente y del vettore casuale
   double r_2;       // Quadrato della distanza dall'origine

   do {

      x_dir = 2.0 * rnd.Rannyu() - 1.0;   // Numero uniforme in [-1,1)
      y_dir = 2.0 * rnd.Rannyu() - 1.0;   // Numero uniforme in [-1,1)

      r_2 = x_dir * x_dir + y_dir * y_dir;   // Distanza quadrata dall'origine

   } while (r_2 >= 1.0 || r_2 == 0.0);
   // Accetto solo punti interni al disco unitario, esclusa l'origine

   double norma = sqrt(r_2);   // Lunghezza del vettore generato

   return y_dir / norma;       // Componente y del versore casuale
}


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
   // PARAMETRI DELL'ESPERIMENTO DI BUFFON
   // =================================================================

   const double lunghezza_ago = 1.0;                    // Lunghezza totale dell'ago
   const double distanza_linee = 2.0 * lunghezza_ago;   // Distanza tra linee consecutive

   const long long numero_lanci = 100000000LL;          // Numero totale di lanci: 10^8
   const int lanci_per_blocco = 100000;                  // Numero di lanci per blocco: 10^5

   const int numero_blocchi =
      numero_lanci / lanci_per_blocco;                  // Numero totale di blocchi: 10^3


   cout << "Numero totale di lanci: " << numero_lanci << endl;
   cout << "Lanci per blocco:       " << lanci_per_blocco << endl;
   cout << "Numero di blocchi:      " << numero_blocchi << endl;


   // =================================================================
   // VETTORI PER IL DATA BLOCKING
   // =================================================================

   vector<double> pi_blocco(numero_blocchi);           // Stima di pi nel singolo blocco
   vector<double> pi_blocco_2(numero_blocchi);         // Quadrato della stima di blocco

   vector<double> pi_progressivo(numero_blocchi);      // Stima progressiva di pi
   vector<double> pi_2_progressivo(numero_blocchi);    // Media progressiva dei quadrati

   vector<double> errore_pi(numero_blocchi);           // Errore statistico progressivo


   // =================================================================
   // SIMULAZIONE DELL'ESPERIMENTO
   // =================================================================

   for (int i = 0; i < numero_blocchi; i++) {

      int numero_hit = 0;   // Numero di intersezioni nel blocco corrente


      for (int j = 0; j < lanci_per_blocco; j++) {

         // -------------------------------------------------------------
         // Posizione casuale del centro dell'ago
         // -------------------------------------------------------------

         double y_centro =
            rnd.Rannyu() * distanza_linee;
         // Il centro dell'ago viene estratto uniformemente in [0,d)


         double distanza_minima =
            min(y_centro, distanza_linee - y_centro);
         // Distanza del centro dell'ago dalla linea piu' vicina


         // -------------------------------------------------------------
         // Orientazione casuale dell'ago
         // -------------------------------------------------------------

         double u_y =
            direzione_y_casuale(rnd);
         // Componente verticale del versore che descrive l'orientazione


         double proiezione_verticale =
            0.5 * lunghezza_ago * fabs(u_y);
         // Proiezione verticale di meta' della lunghezza dell'ago


         // -------------------------------------------------------------
         // Controllo dell'intersezione con una linea
         // -------------------------------------------------------------

         if (proiezione_verticale >= distanza_minima) {

            numero_hit++;   // L'ago interseca una delle linee
         }
      }


      // =================================================================
      // STIMA DI PI NEL BLOCCO
      // =================================================================

      if (numero_hit == 0) {

         cerr << "ERRORE: nessuna intersezione nel blocco "
              << i + 1 << endl;

         return 1;
      }


      pi_blocco[i] =
         (2.0 * lunghezza_ago * lanci_per_blocco) /
         (distanza_linee * numero_hit);
      // pi = 2 L N_thr / (d N_hit)


      pi_blocco_2[i] =
         pi_blocco[i] * pi_blocco[i];
      // Quadrato della stima utilizzato per calcolare l'errore
   }


   // =================================================================
   // MEDIA PROGRESSIVA ED ERRORE STATISTICO
   // =================================================================

   double somma_progressiva = 0.0;      // Somma delle stime di pi dei blocchi
   double somma_2_progressiva = 0.0;    // Somma dei quadrati delle stime


   for (int i = 0; i < numero_blocchi; i++) {

      somma_progressiva += pi_blocco[i];       // Aggiungo la nuova stima
      somma_2_progressiva += pi_blocco_2[i];   // Aggiungo il suo quadrato


      pi_progressivo[i] =
         somma_progressiva / (i + 1);
      // Media delle stime ottenute nei primi i+1 blocchi


      pi_2_progressivo[i] =
         somma_2_progressiva / (i + 1);
      // Media progressiva dei quadrati


      errore_pi[i] =
         errore_progressivo(
            pi_progressivo[i],
            pi_2_progressivo[i],
            i + 1
         );
      // Incertezza statistica sulla stima progressiva
   }


   // =================================================================
   // SALVATAGGIO DEI RISULTATI
   // =================================================================

   ofstream file_output("Pi_Buffon.txt");   // Creo il file di output

   if (!file_output.is_open()) {

      cerr << "ERRORE: impossibile creare Pi_Buffon.txt" << endl;
      return 1;
   }


   file_output << "# Esercizio 01.3 - Esperimento dell'ago di Buffon" << endl;
   file_output << "# Lunghezza ago = " << lunghezza_ago << endl;
   file_output << "# Distanza linee = " << distanza_linee << endl;
   file_output << "# Numero totale di lanci = " << numero_lanci << endl;
   file_output << "# Numero di blocchi = " << numero_blocchi << endl;
   file_output << "# Lanci per blocco = " << lanci_per_blocco << endl;

   file_output << "# Colonna 1: stima di pi nel singolo blocco" << endl;
   file_output << "# Colonna 2: stima progressiva di pi" << endl;
   file_output << "# Colonna 3: errore statistico progressivo" << endl;


   for (int i = 0; i < numero_blocchi; i++) {

      file_output << pi_blocco[i] << " "
                  << pi_progressivo[i] << " "
                  << errore_pi[i] << endl;
      // Scrivo una riga per ogni blocco
   }


   file_output.close();   // Chiudo il file di output


   // =================================================================
   // RISULTATO FINALE A TERMINALE
   // =================================================================

   cout << endl;

   cout << "Stima finale di pi = "
        << pi_progressivo[numero_blocchi - 1]
        << " +/- "
        << errore_pi[numero_blocchi - 1]
        << endl;


   // =================================================================
   // CHIUSURA DEL PROGRAMMA
   // =================================================================

   rnd.SaveSeed();   // Salvo lo stato finale del generatore

   return 0;         // Fine corretta del programma
}