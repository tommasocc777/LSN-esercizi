# Esercizio 10



### Es. 10.1

Compilazione:

mpicxx -std=c++17 -O3 -Wall -Wextra main.cpp random.cpp -o main.exe


Esecuzione:

NPROC=10 bash run_10_1_np.sh


### Es. 10.2

Compilazione:

mpicxx -std=c++17 -O3 -Wall -Wextra main.cpp random.cpp -o main.exe


Esecuzione:

NPROC=10 bash run_trials_np.sh


Analisi:

python analyze_trials_np.py --nproc 10

