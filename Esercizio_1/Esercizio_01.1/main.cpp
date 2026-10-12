/***********************************************************************
 Esercizio 01.1 - Test del generatore di numeri pseudo-casuali

 1) Stima di <r> = 1/2 per r uniforme in [0,1), usando il data blocking.
 2) Stima della varianza sigma^2 = <(r - 1/2)^2> = 1/12.
 3) Test dell'uniformita' tramite il chi-quadro di Pearson.

 File di output:
    AV_prog_Err_prog.txt
    AV_prog_VAR_Err_prog_VAR.txt
    Chi2.txt
************************************************************************/

#include <iostream>     
#include <fstream>      
#include <string>      
#include <vector>       
#include <cmath>      
#include "random.h"     

using namespace std;


// Errore statistico progressivo: sqrt[(<A^2> - <A>^2)/(n_blocchi - 1)]
double errore_progressivo(double media, double media_quadrata, int n_blocchi) {

   if (n_blocchi <= 1) return 0.0;   // Con un solo blocco l'errore non e' stimabile

   double varianza = (media_quadrata - media * media) / (n_blocchi - 1);   // Varianza della media

   if (varianza < 0.0) varianza = 0.0;   // Protezione da errori numerici di arrotondamento

   return sqrt(varianza);   // Restituisce l'incertezza statistica
}


int main() {

   // =================================================================
   // INIZIALIZZAZIONE DEL GENERATORE PSEUDO-CASUALE
   // =================================================================

   Random rnd;                       // Generatore pseudo-casuale
   int seed[4];                      // Quattro interi che definiscono il seed
   int primo_1, primo_2;             // Numeri primi usati dal generatore

   ifstream file_primes("Primes");   // Apro il file contenente i numeri primi

   if (file_primes.is_open()) {
      file_primes >> primo_1 >> primo_2;   // Leggo i due numeri primi
   }
   else {
      cerr << "ERRORE: impossibile aprire il file Primes" << endl;
      return 1;
   }

   file_primes.close();   // Chiudo il file Primes


   ifstream file_seed("seed.in");    // Apro il file contenente il seed
   string proprieta;                 // Stringa usata per cercare RANDOMSEED
   bool seed_trovato = false;        // Controlla se il seed e' stato trovato

   if (file_seed.is_open()) {

      while (file_seed >> proprieta) {   // Leggo il file parola per parola

         if (proprieta == "RANDOMSEED") {

            file_seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];   // Leggo il seed
            rnd.SetRandom(seed, primo_1, primo_2);                   // Inizializzo il generatore
            seed_trovato = true;                                    // Segnalo che il seed e' stato trovato
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
   // ESERCIZIO 01.1.1
   // Stima di <r> = 1/2 mediante data blocking
   // =================================================================

   const int numero_totale = 100000;                         // Numero totale di estrazioni
   const int numero_blocchi = 100;                           // Numero di blocchi
   const int lunghezza_blocco = numero_totale / numero_blocchi;   // Estrazioni per blocco

   vector<double> media_blocco(numero_blocchi);               // Media di r nel singolo blocco
   vector<double> media_blocco_2(numero_blocchi);             // Quadrato della media di blocco
   vector<double> media_progressiva(numero_blocchi);          // Media progressiva di <r>
   vector<double> media_2_progressiva(numero_blocchi);        // Media progressiva dei quadrati
   vector<double> errore_media(numero_blocchi);               // Errore statistico progressivo


   // Calcolo delle medie dei singoli blocchi
   for (int i = 0; i < numero_blocchi; i++) {

      double somma_blocco = 0.0;   // Somma dei valori generati nel blocco corrente

      for (int j = 0; j < lunghezza_blocco; j++) {

         double r = rnd.Rannyu();   // Numero pseudo-casuale uniforme in [0,1)
         somma_blocco += r;         // Aggiorno la somma del blocco
      }

      media_blocco[i] = somma_blocco / lunghezza_blocco;        // Media del blocco i
      media_blocco_2[i] = media_blocco[i] * media_blocco[i];    // Quadrato della media del blocco
   }


   // Calcolo delle medie progressive e degli errori
   double somma_progressiva = 0.0;      // Somma delle medie dei blocchi
   double somma_2_progressiva = 0.0;    // Somma dei quadrati delle medie dei blocchi

   for (int i = 0; i < numero_blocchi; i++) {

      somma_progressiva += media_blocco[i];        // Aggiungo la media del nuovo blocco
      somma_2_progressiva += media_blocco_2[i];    // Aggiungo il quadrato della media

      media_progressiva[i] = somma_progressiva / (i + 1);       // Media sui primi i+1 blocchi
      media_2_progressiva[i] = somma_2_progressiva / (i + 1);   // Media progressiva dei quadrati

      errore_media[i] = errore_progressivo(
         media_progressiva[i],
         media_2_progressiva[i],
         i + 1
      );   // Errore statistico sui primi i+1 blocchi
   }


   // Salvataggio dei risultati della media
   ofstream file_output("AV_prog_Err_prog.txt");   // Creo il file di output

   if (!file_output.is_open()) {
      cerr << "ERRORE: impossibile creare AV_prog_Err_prog.txt" << endl;
      return 1;
   }

   file_output << "# Esercizio 01.1.1 - Stima di <r>" << endl;
   file_output << "# Variabile casuale uniforme r in [0,1)" << endl;
   file_output << "# Valore esatto: <r> = 1/2 = 0.5" << endl;
   file_output << "# Numero totale = " << numero_totale
               << ", blocchi = " << numero_blocchi
               << ", lunghezza blocco = " << lunghezza_blocco << endl;
   file_output << "# Colonna 1: media del singolo blocco" << endl;
   file_output << "# Colonna 2: media progressiva" << endl;
   file_output << "# Colonna 3: errore statistico progressivo" << endl;

   for (int i = 0; i < numero_blocchi; i++) {

      file_output << media_blocco[i] << " "
                  << media_progressiva[i] << " "
                  << errore_media[i] << endl;   // Scrivo i risultati del blocco i
   }

   file_output.close();   // Chiudo il file


   // =================================================================
   // ESERCIZIO 01.1.2
   // Stima di sigma^2 = <(r - 1/2)^2> = 1/12
   // =================================================================

   vector<double> varianza_blocco(numero_blocchi);             // Stima della varianza nel blocco
   vector<double> varianza_blocco_2(numero_blocchi);           // Quadrato della stima di blocco
   vector<double> varianza_progressiva(numero_blocchi);        // Stima progressiva della varianza
   vector<double> varianza_2_progressiva(numero_blocchi);      // Media progressiva dei quadrati
   vector<double> errore_varianza(numero_blocchi);             // Errore statistico progressivo


   // Calcolo di <(r - 1/2)^2> nei singoli blocchi
   for (int i = 0; i < numero_blocchi; i++) {

      double somma_blocco = 0.0;   // Somma di (r - 1/2)^2 nel blocco

      for (int j = 0; j < lunghezza_blocco; j++) {

         double r = rnd.Rannyu();       // Numero uniforme in [0,1)
         double scarto = r - 0.5;       // Scarto rispetto al valore medio teorico

         somma_blocco += scarto * scarto;   // Aggiungo (r - 1/2)^2
      }

      varianza_blocco[i] = somma_blocco / lunghezza_blocco;                    // Stima nel blocco
      varianza_blocco_2[i] = varianza_blocco[i] * varianza_blocco[i];          // Quadrato della stima
   }


   // Calcolo della varianza progressiva e del relativo errore
   somma_progressiva = 0.0;       // Azzero la somma progressiva
   somma_2_progressiva = 0.0;     // Azzero la somma dei quadrati

   for (int i = 0; i < numero_blocchi; i++) {

      somma_progressiva += varianza_blocco[i];         // Aggiungo la nuova stima di blocco
      somma_2_progressiva += varianza_blocco_2[i];     // Aggiungo il suo quadrato

      varianza_progressiva[i] = somma_progressiva / (i + 1);       // Varianza progressiva
      varianza_2_progressiva[i] = somma_2_progressiva / (i + 1);   // Media progressiva dei quadrati

      errore_varianza[i] = errore_progressivo(
         varianza_progressiva[i],
         varianza_2_progressiva[i],
         i + 1
      );   // Errore statistico sulla stima progressiva
   }


   // Salvataggio dei risultati della varianza
   file_output.open("AV_prog_VAR_Err_prog_VAR.txt");   // Apro il secondo file di output

   if (!file_output.is_open()) {
      cerr << "ERRORE: impossibile creare AV_prog_VAR_Err_prog_VAR.txt" << endl;
      return 1;
   }

   file_output << "# Esercizio 01.1.2 - Stima della varianza" << endl;
   file_output << "# Osservabile: (r - 1/2)^2" << endl;
   file_output << "# Valore esatto: sigma^2 = 1/12 = "
               << 1.0 / 12.0 << endl;
   file_output << "# Numero totale = " << numero_totale
               << ", blocchi = " << numero_blocchi
               << ", lunghezza blocco = " << lunghezza_blocco << endl;
   file_output << "# Colonna 1: stima nel singolo blocco" << endl;
   file_output << "# Colonna 2: stima progressiva" << endl;
   file_output << "# Colonna 3: errore statistico progressivo" << endl;

   for (int i = 0; i < numero_blocchi; i++) {

      file_output << varianza_blocco[i] << " "
                  << varianza_progressiva[i] << " "
                  << errore_varianza[i] << endl;   // Scrivo i risultati del blocco i
   }

   file_output.close();   // Chiudo il file


   // =================================================================
   // ESERCIZIO 01.1.3
   // Test del chi-quadro di Pearson sull'uniformita' del generatore
   // =================================================================

   const int numero_bin = 100;               // Numero di intervalli in [0,1)
   const int numero_test = 100;              // Numero di test indipendenti
   const int estrazioni_per_test = 10000;    // Numeri casuali generati per ogni test

   const double valore_atteso =
      static_cast<double>(estrazioni_per_test) / numero_bin;   // Eventi attesi per bin (=100)

   vector<double> chi_quadro(numero_test);    // Valore del chi-quadro per ogni test


   // Eseguo i test del chi-quadro
   for (int i = 0; i < numero_test; i++) {

      vector<int> conteggi(numero_bin, 0);   // Conteggi osservati nei vari bin


      // Genero i numeri casuali e li assegno ai rispettivi bin
      for (int j = 0; j < estrazioni_per_test; j++) {

         double r = rnd.Rannyu();   // Numero uniforme in [0,1)

         int indice_bin =
            static_cast<int>(r * numero_bin);   // Es.: r=0.374 -> bin 37

         if (indice_bin >= numero_bin)
            indice_bin = numero_bin - 1;   // Protezione numerica per un eventuale r=1

         conteggi[indice_bin]++;   // Incremento il numero di eventi nel bin
      }


      double chi_quadro_corrente = 0.0;   // Chi-quadro del test corrente


      // Calcolo chi^2 = somma_k (O_k - E_k)^2 / E_k
      for (int k = 0; k < numero_bin; k++) {

         double osservato = conteggi[k];                       // O_k: eventi osservati nel bin
         double differenza = osservato - valore_atteso;        // O_k - E_k

         chi_quadro_corrente +=
            differenza * differenza / valore_atteso;           // Contributo del bin al chi-quadro
      }

      chi_quadro[i] = chi_quadro_corrente;   // Salvo il risultato del test
   }


   // Salvataggio dei risultati del chi-quadro
   file_output.open("Chi2.txt");   // Apro il file di output

   if (!file_output.is_open()) {
      cerr << "ERRORE: impossibile creare Chi2.txt" << endl;
      return 1;
   }

   file_output << "# Esercizio 01.1.3 - Test del chi-quadro di Pearson" << endl;
   file_output << "# Distribuzione uniforme in [0,1)" << endl;
   file_output << "# Numero di bin = " << numero_bin << endl;
   file_output << "# Numero di test = " << numero_test << endl;
   file_output << "# Estrazioni per test = " << estrazioni_per_test << endl;
   file_output << "# Valore atteso per bin = " << valore_atteso << endl;
   file_output << "# Gradi di liberta' = " << numero_bin - 1 << endl;
   file_output << "# Valore medio atteso di chi^2 = " << numero_bin - 1 << endl;
   file_output << "# Colonna 1: indice del test" << endl;
   file_output << "# Colonna 2: valore del chi-quadro" << endl;

   for (int i = 0; i < numero_test; i++) {

      file_output << i + 1 << " "
                  << chi_quadro[i] << endl;   // Scrivo indice e valore del test
   }

   file_output.close();   // Chiudo il file


   // =================================================================
   // CHIUSURA DEL PROGRAMMA
   // =================================================================

   rnd.SaveSeed();   // Salvo lo stato finale del generatore in seed.out

   return 0;         // Fine del programma
}