ESERCIZIO 10 — NUOVE SIMULAZIONI CON 10 PROCESSI MPI
===================================================

Questi file si usano con IL MAIN.CPP / MAIN.EXE dell'ESERCIZIO 10.2,
che gia' supporta sia --geometry italy sia --geometry both.
Non serve ricompilare il C++ se il main.exe attuale funziona.

1) Copiare nella cartella Esercizio_10.2, insieme a main.exe,
   Primes, seed.in e cap_prov_ita.dat:
       run_10_1_np.sh
       run_trials_np.sh
       analyze_trials_np.py

2) Da MSYS2 MINGW64, per rifare l'ESERCIZIO 10.1 con 10 rank:
       NPROC=10 bash run_10_1_np.sh

   Questo esegue due test (independent e migration ogni 50 generazioni),
   ciascuno su cerchio e quadrato, 1000 generazioni, 500 individui/rank.
   I risultati finiscono in results_10_1_np10/.

3) Per rifare tutti i trial dell'ESERCIZIO 10.2 con 10 rank:
       NPROC=10 bash run_trials_np.sh

   Questo esegue 100 test complessivi: 20 independent e 20 per ciascun
   intervallo di migrazione 10, 20, 50, 100.
   I risultati finiscono in results_10_2_np10/.

4) Per analizzare i risultati 10.2 (da Anaconda Prompt, se necessario):
       py -3 analyze_trials_np.py --nproc 10
   Oppure:
       python analyze_trials_np.py --nproc 10

   I tre CSV si troveranno in results_10_2_np10/.

5) Se in futuro si vogliono utilizzare 11 rank:
       NPROC=11 bash run_10_1_np.sh
       NPROC=11 bash run_trials_np.sh
       py -3 analyze_trials_np.py --nproc 11

   Per 11 rank, un continente e' escluso da ciascun evento di migrazione;
   le coppie vengono rigenerate casualmente al turno successivo.

NOTA IMPORTANTE SULL'HARDWARE:
10 processi MPI non richiedono 10 core fisici, ma se i core disponibili
sono meno di 10, il sovraccarico potrebbe rendere le esecuzioni piu' lente.
Verificare Gestione attivita' > Prestazioni > CPU (core e processori logici).

I risultati dei 4 processi in results_10_2/ NON vengono modificati.
In caso di interruzione lo script riprende dalle prove gia' concluse.
Per rifare una prova con gli stessi parametri, usare FORCE=1.
Se si cambiano parametri comuni (es. generazioni), rinominare prima
la relativa cartella dei risultati, per non mescolare i test.

ATTENZIONE SUL CONFRONTO:
Con 10 rank e 500 individui/rank, si valutano 5000 individui per
continente-generazione aggregata (vs 2000 con 4 rank). Non e' quindi
un confronto a costo computazionale totale costante. Il confronto
independent vs migration A PARITA' di 10 rank rimane invece equo.
