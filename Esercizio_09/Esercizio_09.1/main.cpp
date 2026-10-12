#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <cstdlib>

#include "random.h"

using namespace std;


// ============================================================
// Strutture dati
// ============================================================

// Coordinate cartesiane di una città
struct City {
    double x;
    double y;
};


// Un individuo della popolazione:
// - path contiene l'ordine delle città
// - length contiene la lunghezza totale del percorso
struct Individual {
    vector<int> path;
    double length = 0.0;
};


// Parametri del Genetic Algorithm
struct GAParameters {

    // Numero di città
    int n_cities = 34;

    // Numero di individui nella popolazione
    int population_size = 500;

    // Numero totale di generazioni
    int generations = 1000;

    // Esponente utilizzato nella selezione:
    // j = int(M * r^p)
    double selection_power = 2.0;

    // Probabilità di crossover
    double p_crossover = 0.80;

    // Probabilità delle quattro mutazioni
    double p_pair_swap  = 0.05;
    double p_shift      = 0.05;
    double p_block_swap = 0.05;
    double p_inversion  = 0.05;

    // Numero di individui migliori conservati
    // direttamente nella generazione successiva
    int n_elite = 0;

    // Generazioni nelle quali viene salvato
    // il miglior percorso presente nella popolazione
    vector<int> snapshot_generations = {
        0, 10, 20, 50, 100, 200, 500, 1000
    };
};


// ============================================================
// Generatore di numeri interi casuali
// ============================================================

// Estrae un intero uniforme nell'intervallo [min, max]
int RandomInt(Random& rnd, int min, int max) {

    int value = min + static_cast<int>(
        rnd.Rannyu() * (max - min + 1)
    );

    // Protezione numerica nel caso limite
    if (value > max)
        value = max;

    return value;
}


// ============================================================
// Inizializzazione del generatore random
// ============================================================

void InitializeRandom(Random& rnd) {

    int seed[4];
    int p1, p2;

    // Lettura dei numeri primi
    ifstream primes("Primes");

    if (!primes.is_open()) {
        cerr << "ERRORE: impossibile aprire il file Primes." << endl;
        exit(1);
    }

    primes >> p1 >> p2;
    primes.close();


    // Lettura del seed
    ifstream input("seed.in");

    if (!input.is_open()) {
        cerr << "ERRORE: impossibile aprire il file seed.in." << endl;
        exit(1);
    }

    string property;
    bool seed_found = false;

    while (input >> property) {

        if (property == "RANDOMSEED") {

            input >> seed[0]
                  >> seed[1]
                  >> seed[2]
                  >> seed[3];

            rnd.SetRandom(seed, p1, p2);

            seed_found = true;
            break;
        }
    }

    input.close();

    if (!seed_found) {
        cerr << "ERRORE: RANDOMSEED non trovato in seed.in." << endl;
        exit(1);
    }
}


// ============================================================
// Generazione delle città
// ============================================================

// ------------------------------------------------------------
// Città distribuite sulla circonferenza unitaria
// ------------------------------------------------------------

vector<City> GenerateCitiesCircle(Random& rnd, int n_cities) {

    vector<City> cities(n_cities);

    const double PI = acos(-1.0);

    for (int i = 0; i < n_cities; ++i) {

        // Angolo uniforme in [0, 2*pi)
        double theta = rnd.Rannyu(0.0, 2.0 * PI);

        cities[i].x = cos(theta);
        cities[i].y = sin(theta);
    }

    return cities;
}


// ------------------------------------------------------------
// Città distribuite uniformemente nel quadrato unitario
// ------------------------------------------------------------

vector<City> GenerateCitiesSquare(Random& rnd, int n_cities) {

    vector<City> cities(n_cities);

    for (int i = 0; i < n_cities; ++i) {

        cities[i].x = rnd.Rannyu();
        cities[i].y = rnd.Rannyu();
    }

    return cities;
}


// ============================================================
// Creazione di un cromosoma casuale
// ============================================================

vector<int> RandomPath(Random& rnd, int n_cities) {

    vector<int> path(n_cities);

    // path = [1, 2, 3, ..., N]
    iota(path.begin(), path.end(), 1);


    // Fisher-Yates shuffle.
    //
    // La posizione 0 non viene modificata:
    // la città 1 rimane quindi sempre la prima città.
    for (int i = n_cities - 1; i >= 2; --i) {

        int j = RandomInt(rnd, 1, i);

        swap(path[i], path[j]);
    }

    return path;
}


// ============================================================
// Controllo della validità di un percorso
// ============================================================

bool CheckPath(const vector<int>& path, int n_cities) {

    // Il cromosoma deve contenere esattamente N città
    if (static_cast<int>(path.size()) != n_cities)
        return false;


    // La prima città deve essere sempre la città 1
    if (path[0] != 1)
        return false;


    // Contatore delle occorrenze delle città
    vector<int> count(n_cities + 1, 0);

    for (int city : path) {

        // Ogni indice deve appartenere a [1, N]
        if (city < 1 || city > n_cities)
            return false;

        count[city]++;
    }


    // Ogni città deve comparire esattamente una volta
    for (int city = 1; city <= n_cities; ++city) {

        if (count[city] != 1)
            return false;
    }

    return true;
}


// ============================================================
// Calcolo della funzione costo L^(1)
// ============================================================

double PathLength(
    const vector<int>& path,
    const vector<City>& cities
) {

    double length = 0.0;

    int n_cities = static_cast<int>(path.size());


    for (int i = 0; i < n_cities; ++i) {

        // Città corrente
        int current = path[i] - 1;

        // Città successiva.
        // Il modulo chiude automaticamente il percorso.
        int next = path[(i + 1) % n_cities] - 1;


        double dx = cities[current].x - cities[next].x;
        double dy = cities[current].y - cities[next].y;


        // Distanza euclidea
        length += sqrt(dx * dx + dy * dy);
    }

    return length;
}


// ============================================================
// Valutazione e ordinamento della popolazione
// ============================================================

void EvaluatePopulation(
    vector<Individual>& population,
    const vector<City>& cities
) {

    // Calcolo della funzione costo di ogni individuo
    for (Individual& individual : population) {

        individual.length =
            PathLength(individual.path, cities);
    }


    // Ordinamento crescente:
    // il primo individuo è quello con L più piccolo
    sort(
        population.begin(),
        population.end(),
        [](const Individual& a, const Individual& b) {
            return a.length < b.length;
        }
    );
}


// ============================================================
// Media di L sulla migliore metà della popolazione
// ============================================================

double BestHalfMean(const vector<Individual>& population) {

    int half = static_cast<int>(population.size()) / 2;

    double sum = 0.0;

    for (int i = 0; i < half; ++i)
        sum += population[i].length;

    return sum / static_cast<double>(half);
}


// ============================================================
// Selezione
// ============================================================

// La popolazione deve essere già ordinata secondo L.
//
// j = int(M * r^p)
//
// Con p > 1 vengono favoriti gli indici piccoli,
// cioè gli individui con funzione costo minore.
int SelectIndividual(
    Random& rnd,
    int population_size,
    double selection_power
) {

    double r = rnd.Rannyu();

    int index = static_cast<int>(
        population_size * pow(r, selection_power)
    );


    // Protezione numerica
    if (index >= population_size)
        index = population_size - 1;

    return index;
}


// ============================================================
// MUTAZIONE 1:
// scambio di due città
// ============================================================

void MutationPairSwap(
    Individual& individual,
    Random& rnd
) {

    int n = static_cast<int>(individual.path.size());


    // La posizione 0 non viene mai modificata
    int i = RandomInt(rnd, 1, n - 1);
    int j = RandomInt(rnd, 1, n - 1);


    // Vengono scelte due posizioni differenti
    while (j == i)
        j = RandomInt(rnd, 1, n - 1);


    swap(
        individual.path[i],
        individual.path[j]
    );
}


// ============================================================
// MUTAZIONE 2:
// shift di un blocco
// ============================================================

void MutationShift(
    Individual& individual,
    Random& rnd
) {

    int n = static_cast<int>(individual.path.size());

    // Numero di città modificabili
    int movable = n - 1;


    // Lunghezza del blocco da traslare
    int m = RandomInt(rnd, 1, movable - 1);


    // Numero di posizioni dello shift
    int shift = RandomInt(rnd, 1, movable - m);


    // Lunghezza complessiva della regione coinvolta
    int total = m + shift;


    // Posizione iniziale della regione
    int start = RandomInt(rnd, 1, n - total);


    // Esempio:
    //
    // [1, 2, 3, 4, 5]
    //
    // m = 2, shift = 2
    //
    // [1, 4, 5, 2, 3]
    rotate(
        individual.path.begin() + start,
        individual.path.begin() + start + m,
        individual.path.begin() + start + total
    );
}


// ============================================================
// MUTAZIONE 3:
// scambio di due blocchi della stessa lunghezza
// ============================================================

void MutationBlockSwap(
    Individual& individual,
    Random& rnd
) {

    int n = static_cast<int>(individual.path.size());

    int movable = n - 1;


    // La lunghezza dei blocchi non può superare
    // metà delle città modificabili
    int m = RandomInt(rnd, 1, movable / 2);


    // Posizione iniziale del primo blocco
    int first_start = RandomInt(rnd, 1, n - m);


    // Possibili posizioni del secondo blocco
    // che non si sovrappongono al primo
    vector<int> possible_starts;


    for (int second_start = 1;
         second_start <= n - m;
         ++second_start) {

        bool non_overlapping =
            (second_start + m - 1 < first_start) ||
            (first_start + m - 1 < second_start);

        if (non_overlapping)
            possible_starts.push_back(second_start);
    }


    // In caso non esista un secondo blocco valido
    // la mutazione non viene applicata
    if (possible_starts.empty())
        return;


    int chosen =
        RandomInt(
            rnd,
            0,
            static_cast<int>(possible_starts.size()) - 1
        );

    int second_start = possible_starts[chosen];


    // Scambio elemento per elemento
    for (int k = 0; k < m; ++k) {

        swap(
            individual.path[first_start + k],
            individual.path[second_start + k]
        );
    }
}


// ============================================================
// MUTAZIONE 4:
// inversione di un blocco
// ============================================================

void MutationInversion(
    Individual& individual,
    Random& rnd
) {

    int n = static_cast<int>(individual.path.size());

    int movable = n - 1;


    // Lunghezza della sequenza da invertire
    int m = RandomInt(rnd, 2, movable);


    // Posizione iniziale
    int start = RandomInt(rnd, 1, n - m);


    reverse(
        individual.path.begin() + start,
        individual.path.begin() + start + m
    );
}


// ============================================================
// Applicazione delle mutazioni
// ============================================================

void ApplyMutations(
    Individual& individual,
    Random& rnd,
    const GAParameters& parameters
) {

    // Le quattro mutazioni vengono testate
    // indipendentemente.

    if (rnd.Rannyu() < parameters.p_pair_swap)
        MutationPairSwap(individual, rnd);


    if (rnd.Rannyu() < parameters.p_shift)
        MutationShift(individual, rnd);


    if (rnd.Rannyu() < parameters.p_block_swap)
        MutationBlockSwap(individual, rnd);


    if (rnd.Rannyu() < parameters.p_inversion)
        MutationInversion(individual, rnd);
}


// ============================================================
// Crossover
// ============================================================

pair<Individual, Individual> Crossover(
    const Individual& parent1,
    const Individual& parent2,
    Random& rnd
) {

    int n = static_cast<int>(parent1.path.size());


    // Punto di taglio.
    //
    // Le prime "cut" città vengono conservate.
    int cut = RandomInt(rnd, 1, n - 1);


    Individual child1;
    Individual child2;


    child1.path.assign(n, -1);
    child2.path.assign(n, -1);


    // ========================================================
    // Costruzione del primo figlio
    // ========================================================

    vector<bool> used1(n + 1, false);


    // Prima parte ereditata dal genitore 1
    for (int i = 0; i < cut; ++i) {

        child1.path[i] = parent1.path[i];

        used1[parent1.path[i]] = true;
    }


    // Le posizioni rimanenti vengono completate
    // seguendo l'ordine del genitore 2
    int position = cut;

    for (int city : parent2.path) {

        if (!used1[city]) {

            child1.path[position] = city;

            used1[city] = true;

            position++;
        }
    }


    // ========================================================
    // Costruzione del secondo figlio
    // ========================================================

    vector<bool> used2(n + 1, false);


    // Prima parte ereditata dal genitore 2
    for (int i = 0; i < cut; ++i) {

        child2.path[i] = parent2.path[i];

        used2[parent2.path[i]] = true;
    }


    // Le posizioni rimanenti vengono completate
    // seguendo l'ordine del genitore 1
    position = cut;

    for (int city : parent1.path) {

        if (!used2[city]) {

            child2.path[position] = city;

            used2[city] = true;

            position++;
        }
    }


    return {child1, child2};
}


// ============================================================
// Creazione della popolazione iniziale
// ============================================================

vector<Individual> CreateInitialPopulation(
    Random& rnd,
    const GAParameters& parameters
) {

    vector<Individual> population;

    population.reserve(parameters.population_size);


    for (int i = 0;
         i < parameters.population_size;
         ++i) {

        Individual individual;

        individual.path =
            RandomPath(rnd, parameters.n_cities);


        // Controllo di sicurezza
        if (!CheckPath(
                individual.path,
                parameters.n_cities
            )) {

            cerr
                << "ERRORE: cromosoma iniziale non valido."
                << endl;

            exit(1);
        }


        population.push_back(individual);
    }

    return population;
}


// ============================================================
// Creazione della generazione successiva
// ============================================================

vector<Individual> NewGeneration(
    const vector<Individual>& population,
    Random& rnd,
    const GAParameters& parameters
) {

    vector<Individual> new_population;

    new_population.reserve(parameters.population_size);


    // ========================================================
    // Elitismo
    // ========================================================

    for (int i = 0;
         i < parameters.n_elite &&
         i < parameters.population_size;
         ++i) {

        new_population.push_back(population[i]);
    }


    // ========================================================
    // Generazione dei nuovi individui
    // ========================================================

    while (
        static_cast<int>(new_population.size())
        < parameters.population_size
    ) {

        // ----------------------------------------------------
        // Selezione di due genitori
        // ----------------------------------------------------

        int index1 =
            SelectIndividual(
                rnd,
                parameters.population_size,
                parameters.selection_power
            );

        int index2 =
            SelectIndividual(
                rnd,
                parameters.population_size,
                parameters.selection_power
            );


        const Individual& parent1 = population[index1];
        const Individual& parent2 = population[index2];


        Individual child1;
        Individual child2;


        // ----------------------------------------------------
        // Crossover
        // ----------------------------------------------------

        if (rnd.Rannyu() < parameters.p_crossover) {

            auto children =
                Crossover(parent1, parent2, rnd);

            child1 = children.first;
            child2 = children.second;
        }

        else {

            // Nessun crossover:
            // i figli sono copie dei genitori
            child1 = parent1;
            child2 = parent2;
        }


        // ----------------------------------------------------
        // Mutazioni
        // ----------------------------------------------------

        ApplyMutations(
            child1,
            rnd,
            parameters
        );

        ApplyMutations(
            child2,
            rnd,
            parameters
        );


        // ----------------------------------------------------
        // Controllo della validità
        // ----------------------------------------------------

        if (!CheckPath(
                child1.path,
                parameters.n_cities
            )) {

            cerr
                << "ERRORE: primo figlio non valido."
                << endl;

            exit(1);
        }


        if (!CheckPath(
                child2.path,
                parameters.n_cities
            )) {

            cerr
                << "ERRORE: secondo figlio non valido."
                << endl;

            exit(1);
        }


        // ----------------------------------------------------
        // Inserimento dei figli
        // ----------------------------------------------------

        new_population.push_back(child1);


        // Se rimane ancora spazio viene inserito
        // anche il secondo figlio
        if (
            static_cast<int>(new_population.size())
            < parameters.population_size
        ) {

            new_population.push_back(child2);
        }
    }


    return new_population;
}


// ============================================================
// Salvataggio delle città
// ============================================================

void SaveCities(
    const vector<City>& cities,
    const string& filename
) {

    ofstream output(filename);

    if (!output.is_open()) {

        cerr
            << "ERRORE: impossibile creare "
            << filename
            << endl;

        exit(1);
    }


    output << "# Coordinate delle citta'\n";
    output << "#\n";
    output << "# Colonna 1: indice della citta'\n";
    output << "# Colonna 2: coordinata x\n";
    output << "# Colonna 3: coordinata y\n";
    output << "#\n";
    output << "# city x y\n";


    output << scientific << setprecision(12);


    for (size_t i = 0; i < cities.size(); ++i) {

        output
            << i + 1 << " "
            << cities[i].x << " "
            << cities[i].y
            << "\n";
    }


    output.close();
}


// ============================================================
// Salvataggio di un percorso
// ============================================================

void SaveBestPath(
    const Individual& best,
    const vector<City>& cities,
    const string& filename
) {

    ofstream output(filename);

    if (!output.is_open()) {

        cerr
            << "ERRORE: impossibile creare "
            << filename
            << endl;

        exit(1);
    }


    output
        << "# Miglior percorso trovato dal Genetic Algorithm\n";

    output << "#\n";

    output
        << "# L_best = "
        << scientific
        << setprecision(12)
        << best.length
        << "\n";

    output << "#\n";
    output << "# Colonna 1: ordine di visita\n";
    output << "# Colonna 2: indice della citta'\n";
    output << "# Colonna 3: coordinata x\n";
    output << "# Colonna 4: coordinata y\n";
    output << "#\n";
    output << "# L'ultima riga ripete la prima citta'\n";
    output << "# per chiudere graficamente il percorso.\n";
    output << "#\n";
    output << "# order city x y\n";


    output << scientific << setprecision(12);


    // Percorso completo
    for (size_t i = 0; i < best.path.size(); ++i) {

        int city = best.path[i];

        output
            << i + 1 << " "
            << city << " "
            << cities[city - 1].x << " "
            << cities[city - 1].y
            << "\n";
    }


    // Ripetizione della prima città
    int first_city = best.path[0];

    output
        << best.path.size() + 1 << " "
        << first_city << " "
        << cities[first_city - 1].x << " "
        << cities[first_city - 1].y
        << "\n";


    output.close();
}


// ============================================================
// Controllo se una generazione deve essere salvata
// ============================================================

bool IsSnapshotGeneration(
    int generation,
    const vector<int>& snapshot_generations
) {

    return find(
        snapshot_generations.begin(),
        snapshot_generations.end(),
        generation
    ) != snapshot_generations.end();
}


// ============================================================
// Salvataggio README
// ============================================================

void WriteReadme(
    const GAParameters& parameters
) {

    ofstream output("README_09.txt");

    if (!output.is_open()) {

        cerr
            << "ERRORE: impossibile creare README_09.txt."
            << endl;

        exit(1);
    }


    output << "ESERCIZIO 09.1 - GENETIC ALGORITHM PER IL TSP\n";
    output << "==============================================\n\n";

    output << "PARAMETRI\n";
    output << "----------\n";
    output << "Numero citta': "
           << parameters.n_cities << "\n";

    output << "Popolazione: "
           << parameters.population_size << "\n";

    output << "Generazioni: "
           << parameters.generations << "\n";

    output << "Esponente selezione p: "
           << parameters.selection_power << "\n";

    output << "Probabilita' crossover: "
           << parameters.p_crossover << "\n";

    output << "Probabilita' pair swap: "
           << parameters.p_pair_swap << "\n";

    output << "Probabilita' shift: "
           << parameters.p_shift << "\n";

    output << "Probabilita' block swap: "
           << parameters.p_block_swap << "\n";

    output << "Probabilita' inversione: "
           << parameters.p_inversion << "\n";

    output << "Numero elite: "
           << parameters.n_elite << "\n\n";


    output << "RAPPRESENTAZIONE\n";
    output << "----------------\n";

    output
        << "Ogni cromosoma e' una permutazione delle citta'.\n";

    output
        << "La citta' 1 viene mantenuta sempre nella prima posizione.\n\n";


    output << "FUNZIONE COSTO\n";
    output << "--------------\n";

    output
        << "L^(1) e' la somma delle distanze euclidee tra citta' consecutive,\n";

    output
        << "inclusa la distanza tra l'ultima citta' e la prima.\n\n";


    output << "SNAPSHOT DEL MIGLIOR PERCORSO\n";
    output << "-----------------------------\n";

    output
        << "Il miglior percorso viene salvato alle generazioni:\n";

    for (int generation :
         parameters.snapshot_generations) {

        output << generation << " ";
    }

    output << "\n\n";


    output << "FILE PRINCIPALI\n";
    output << "---------------\n";

    output << "circle_cities.txt\n";
    output << "circle_convergence.txt\n";
    output << "circle_best_path.txt\n\n";

    output << "square_cities.txt\n";
    output << "square_convergence.txt\n";
    output << "square_best_path.txt\n\n";


    output << "FILE SNAPSHOT\n";
    output << "-------------\n";

    output
        << "circle_best_path_gen_X.txt\n";

    output
        << "square_best_path_gen_X.txt\n";

    output
        << "dove X indica il numero della generazione.\n";


    output.close();
}


// ============================================================
// Esecuzione del Genetic Algorithm
// ============================================================

void RunGeneticAlgorithm(
    Random& rnd,
    const vector<City>& cities,
    const GAParameters& parameters,
    const string& name
) {

    cout << "\n============================================\n";
    cout << "Simulazione: " << name << endl;
    cout << "============================================\n";


    // ========================================================
    // Popolazione iniziale
    // ========================================================

    vector<Individual> population =
        CreateInitialPopulation(
            rnd,
            parameters
        );


    // Calcolo delle funzioni costo
    // e ordinamento della popolazione
    EvaluatePopulation(
        population,
        cities
    );


    // ========================================================
    // File di convergenza
    // ========================================================

    string convergence_filename =
        name + "_convergence.txt";

    ofstream convergence(
        convergence_filename
    );


    if (!convergence.is_open()) {

        cerr
            << "ERRORE: impossibile creare "
            << convergence_filename
            << endl;

        exit(1);
    }


    convergence
        << "# Convergenza del Genetic Algorithm\n";

    convergence << "#\n";

    convergence
        << "# Colonna 1: generazione\n";

    convergence
        << "# Colonna 2: L^(1) del miglior percorso\n";

    convergence
        << "# Colonna 3: <L^(1)> della migliore meta' della popolazione\n";

    convergence << "#\n";
    convergence
        << "# generation L_best L_best_half_mean\n";


    convergence
        << scientific
        << setprecision(12);


    // ========================================================
    // Generazione 0
    // ========================================================

    convergence
        << 0 << " "
        << population[0].length << " "
        << BestHalfMean(population)
        << "\n";


    // Salvataggio del miglior percorso iniziale
    if (
        IsSnapshotGeneration(
            0,
            parameters.snapshot_generations
        )
    ) {

        string filename =
            name
            + "_best_path_gen_0.txt";

        SaveBestPath(
            population[0],
            cities,
            filename
        );
    }


    // ========================================================
    // Evoluzione
    // ========================================================

    for (
        int generation = 1;
        generation <= parameters.generations;
        ++generation
    ) {

        // Creazione della nuova popolazione
        population =
            NewGeneration(
                population,
                rnd,
                parameters
            );


        // Valutazione e ordinamento
        EvaluatePopulation(
            population,
            cities
        );


        // Valori di convergenza
        double best_length =
            population[0].length;

        double best_half_mean =
            BestHalfMean(population);


        convergence
            << generation << " "
            << best_length << " "
            << best_half_mean
            << "\n";


        // ====================================================
        // Salvataggio degli snapshot
        // ====================================================

        if (
            IsSnapshotGeneration(
                generation,
                parameters.snapshot_generations
            )
        ) {

            string filename =
                name
                + "_best_path_gen_"
                + to_string(generation)
                + ".txt";

            SaveBestPath(
                population[0],
                cities,
                filename
            );
        }


        // Informazioni a terminale ogni 100 generazioni
        if (generation % 100 == 0) {

            cout
                << "Generazione "
                << setw(4)
                << generation
                << "  L_best = "
                << best_length
                << endl;
        }
    }


    convergence.close();


    // ========================================================
    // Salvataggio del miglior percorso finale
    // ========================================================

    string best_path_filename =
        name + "_best_path.txt";

    SaveBestPath(
        population[0],
        cities,
        best_path_filename
    );


    // ========================================================
    // Risultato finale a terminale
    // ========================================================

    cout << "\nMiglior percorso finale (" << name << "):\n";

    cout << "L_best = "
         << population[0].length
         << "\n\n";

    cout << "Percorso:\n";

    for (int city : population[0].path)
        cout << city << " -> ";

    // Ritorno alla città iniziale
    cout << population[0].path[0] << endl;
}


// ============================================================
// MAIN
// ============================================================

int main() {

    // ========================================================
    // Inizializzazione del generatore random
    // ========================================================

    Random rnd;

    InitializeRandom(rnd);


    // ========================================================
    // Parametri del Genetic Algorithm
    // ========================================================

    GAParameters parameters;


    // ========================================================
    // Generazione delle due configurazioni geometriche
    //
    // Entrambi gli insiemi di città vengono generati PRIMA
    // delle simulazioni genetiche.
    //
    // In questo modo la geometria del quadrato non dipende
    // dal numero di operazioni casuali effettuate durante
    // il GA sulla circonferenza.
    // ========================================================

    vector<City> circle_cities =
        GenerateCitiesCircle(
            rnd,
            parameters.n_cities
        );

    vector<City> square_cities =
        GenerateCitiesSquare(
            rnd,
            parameters.n_cities
        );


    // ========================================================
    // Salvataggio delle coordinate
    // ========================================================

    SaveCities(
        circle_cities,
        "circle_cities.txt"
    );

    SaveCities(
        square_cities,
        "square_cities.txt"
    );


    // ========================================================
    // README
    // ========================================================

    WriteReadme(parameters);


    // ========================================================
    // Genetic Algorithm: circonferenza
    // ========================================================

    RunGeneticAlgorithm(
    rnd,
    circle_cities,
    parameters,
    "circle_no_elite"
    );


    // ========================================================
    // Genetic Algorithm: quadrato
    // ========================================================

    RunGeneticAlgorithm(
    rnd,
    square_cities,
    parameters,
    "square_no_elite"
    );


    // ========================================================
    // Salvataggio del seed finale
    // ========================================================

    rnd.SaveSeed();


    cout << "\n============================================\n";
    cout << "Simulazioni completate." << endl;
    cout << "============================================\n";


    return 0;
}