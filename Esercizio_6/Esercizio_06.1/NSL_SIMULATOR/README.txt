ESERCIZIO 06.1 - ISING 1D
===========================

Eseguire sempre il programma dalla cartella SOURCE.

Compilazione:
  make clean
  make

Esecuzione:
  ./simulator.exe
  oppure: make run

SCELTE IN INPUT/input.dat
-------------------------
RUN_MODE EQ_ONLY
  Esegue una sola equilibrazione e si ferma.
  Output: equilibration_energy.dat e equilibration_blocks.dat.

RUN_MODE SINGLE
  Esegue equilibrazione + produzione.
  Usa properties.dat per decidere le osservabili da misurare.

RUN_MODE ALL_EQ
  Esegue automaticamente le 12 prove di equilibrazione:
  T = 0.5, 1.0, 2.0; Metropolis/Gibbs; ordered/random.
  I file vengono salvati con i nomi gia' attesi dal notebook Es_06.ipynb.

ALGORITMO
  SIMULATION_TYPE 2 1.0 0.00  -> Metropolis, J=1, h=0
  SIMULATION_TYPE 3 1.0 0.00  -> Gibbs,      J=1, h=0

CONFIGURAZIONE INIZIALE
  RESTART 0 + START_STATE ORDERED -> tutti gli spin +1
  RESTART 0 + START_STATE RANDOM  -> spin casuali +/-1
  RESTART 1                        -> legge INPUT/CONFIG/config.spin
                                      (START_STATE viene ignorato)

TEMPERATURA
  TEMP 1.0 -> T=1.0
  Cambiare solo questo valore per una singola run.

EQUILIBRAZIONE
  NEQBLOCKS 2000
  NEQSTEPS  10
  Totale = 20000 MC steps.
  equilibration_energy.dat salva U/N a ogni MC step, incluso step 0.
  equilibration_blocks.dat salva la media locale ogni 10 MC steps.

PRODUZIONE
  NBLOCKS 100
  NSTEPS  10000
  Totale = 10^6 MC steps.

PROPERTIES
  h=0: TOTAL_ENERGY, SPECIFIC_HEAT, SUSCEPTIBILITY
  h=0.02: MAGNETIZATION

I template pronti sono in INPUT/TEMPLATES.

AGGIORNAMENTO EQUILIBRAZIONE
---------------------------
Durante l'equilibrazione vengono ora salvati, a ogni MC step,
U/N e M/N nel file equilibration_energy*.dat.
Nel file equilibration_blocks*.dat vengono salvate anche le medie locali
BLOCK_U/N e BLOCK_M/N su NEQSTEPS step.
I dati istantanei permettono di costruire nel notebook anche C/N e chi/N
dalle fluttuazioni, senza modificare il data blocking di produzione.
