ESERCIZIO 09.1 - GENETIC ALGORITHM PER IL TSP
==============================================

PARAMETRI
----------
Numero citta': 34
Popolazione: 500
Generazioni: 1000
Esponente selezione p: 2
Probabilita' crossover: 0.8
Probabilita' pair swap: 0.05
Probabilita' shift: 0.05
Probabilita' block swap: 0.05
Probabilita' inversione: 0.05
Numero elite: 0

RAPPRESENTAZIONE
----------------
Ogni cromosoma e' una permutazione delle citta'.
La citta' 1 viene mantenuta sempre nella prima posizione.

FUNZIONE COSTO
--------------
L^(1) e' la somma delle distanze euclidee tra citta' consecutive,
inclusa la distanza tra l'ultima citta' e la prima.

SNAPSHOT DEL MIGLIOR PERCORSO
-----------------------------
Il miglior percorso viene salvato alle generazioni:
0 10 20 50 100 200 500 1000 

FILE PRINCIPALI
---------------
circle_cities.txt
circle_convergence.txt
circle_best_path.txt

square_cities.txt
square_convergence.txt
square_best_path.txt

FILE SNAPSHOT
-------------
circle_best_path_gen_X.txt
square_best_path_gen_X.txt
dove X indica il numero della generazione.
