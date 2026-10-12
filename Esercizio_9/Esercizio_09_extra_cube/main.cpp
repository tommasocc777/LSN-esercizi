#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <utility>
#include <cstdlib>

#include "random.h"

using namespace std;


// ============================================================
// Strutture dati
// ============================================================

// Una città viene sempre rappresentata tramite coordinate 3D.
// Per circonferenza e quadrato viene semplicemente posto z = 0.
struct City {
    double x;
    double y;
    double z;
};


// Individuo della popolazione
struct Individual {
    vector<int> path;
    double length = 0.0;
};


// Geometrie considerate
enum class Geometry {
    Circle,
    Square,
    Sphere,
    Cube
};


// Parametri del Genetic Algorithm
struct GAParameters {

    int n_cities = 100;
    int population_size = 500;
    int generations = 2000;

    // Selezione:
    // j = int(M * r^p)
    double selection_power = 2.0;

    // Probabilità di crossover
    double p_crossover = 0.80;

    // Probabilità delle quattro mutazioni
    double p_pair_swap  = 0.05;
    double p_shift      = 0.05;
    double p_block_swap = 0.05;
    double p_inversion  = 0.05;

    // Elitismo
    int n_elite = 1;

    // Generazioni in cui viene salvato il miglior percorso
    vector<int> snapshot_generations = {
        0, 10, 20, 50, 100, 200, 500, 1000
    };
};


// ============================================================
// Generazione di numeri interi casuali
// ============================================================

int RandomInt(Random& rnd, int min, int max) {

    int value =
        min + static_cast<int>(
            rnd.Rannyu() * (max - min + 1)
        );

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

    ifstream primes("Primes");

    if (!primes.is_open()) {
        cerr << "ERRORE: impossibile aprire Primes." << endl;
        exit(1);
    }

    primes >> p1 >> p2;
    primes.close();


    ifstream input("seed.in");

    if (!input.is_open()) {
        cerr << "ERRORE: impossibile aprire seed.in." << endl;
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
        cerr << "ERRORE: RANDOMSEED non trovato." << endl;
        exit(1);
    }
}


// ============================================================
// Generazione delle città
// ============================================================

// ------------------------------------------------------------
// Circonferenza unitaria
// ------------------------------------------------------------

vector<City> GenerateCitiesCircle(
    Random& rnd,
    int n_cities
) {

    vector<City> cities(n_cities);

    const double PI = acos(-1.0);

    for (int i = 0; i < n_cities; ++i) {

        double theta =
            rnd.Rannyu(0.0, 2.0 * PI);

        cities[i].x = cos(theta);
        cities[i].y = sin(theta);
        cities[i].z = 0.0;
    }

    return cities;
}


// ------------------------------------------------------------
// Quadrato unitario
// ------------------------------------------------------------

vector<City> GenerateCitiesSquare(
    Random& rnd,
    int n_cities
) {

    vector<City> cities(n_cities);

    for (int i = 0; i < n_cities; ++i) {

        cities[i].x = rnd.Rannyu();
        cities[i].y = rnd.Rannyu();
        cities[i].z = 0.0;
    }

    return cities;
}


// ------------------------------------------------------------
// Sfera unitaria
// ------------------------------------------------------------

vector<City> GenerateCitiesSphere(
    Random& rnd,
    int n_cities
) {

    vector<City> cities(n_cities);

    const double PI = acos(-1.0);

    for (int i = 0; i < n_cities; ++i) {

        // Distribuzione uniforme sulla superficie della sfera
        double z = rnd.Rannyu(-1.0, 1.0);
        double phi = rnd.Rannyu(0.0, 2.0 * PI);

        double rho =
            sqrt(1.0 - z * z);

        cities[i].x =
            rho * cos(phi);

        cities[i].y =
            rho * sin(phi);

        cities[i].z = z;
    }

    return cities;
}


// ------------------------------------------------------------
// Cubo unitario
//
// Le città sono distribuite uniformemente nel volume:
// x, y, z appartenenti a [0,1)
// ------------------------------------------------------------

vector<City> GenerateCitiesCube(
    Random& rnd,
    int n_cities
) {

    vector<City> cities(n_cities);

    for (int i = 0; i < n_cities; ++i) {

        cities[i].x = rnd.Rannyu();
        cities[i].y = rnd.Rannyu();
        cities[i].z = rnd.Rannyu();
    }

    return cities;
}


// ============================================================
// Generazione di un cromosoma casuale
// ============================================================

vector<int> RandomPath(
    Random& rnd,
    int n_cities
) {

    vector<int> path(n_cities);

    // [1, 2, 3, ..., N]
    iota(path.begin(), path.end(), 1);


    // Fisher-Yates.
    // La città 1 resta sempre nella prima posizione.
    for (int i = n_cities - 1; i >= 2; --i) {

        int j =
            RandomInt(rnd, 1, i);

        swap(path[i], path[j]);
    }

    return path;
}


// ============================================================
// Controllo della validità del percorso
// ============================================================

bool CheckPath(
    const vector<int>& path,
    int n_cities
) {

    if (static_cast<int>(path.size()) != n_cities)
        return false;


    // La città 1 deve restare nella prima posizione
    if (path[0] != 1)
        return false;


    vector<int> count(n_cities + 1, 0);


    for (int city : path) {

        if (city < 1 || city > n_cities)
            return false;

        count[city]++;
    }


    // Ogni città deve comparire una sola volta
    for (int city = 1;
         city <= n_cities;
         ++city) {

        if (count[city] != 1)
            return false;
    }

    return true;
}


// ============================================================
// Distanze
// ============================================================

// ------------------------------------------------------------
// Distanza euclidea 3D
//
// Viene utilizzata per:
// - circonferenza
// - quadrato
// - cubo
// ------------------------------------------------------------

double EuclideanDistance(
    const City& a,
    const City& b
) {

    double dx = a.x - b.x;
    double dy = a.y - b.y;
    double dz = a.z - b.z;

    return sqrt(
        dx * dx +
        dy * dy +
        dz * dz
    );
}


// ------------------------------------------------------------
// Distanza geodetica sulla sfera unitaria
//
// d = arccos(r_i dot r_j)
// ------------------------------------------------------------

double GeodesicDistance(
    const City& a,
    const City& b
) {

    double dot =
        a.x * b.x +
        a.y * b.y +
        a.z * b.z;


    // Correzione di eventuali errori numerici
    dot =
        max(-1.0, min(1.0, dot));


    return acos(dot);
}


// ============================================================
// Funzione costo
// ============================================================

double PathLength(
    const vector<int>& path,
    const vector<City>& cities,
    Geometry geometry
) {

    double length = 0.0;

    int n =
        static_cast<int>(path.size());


    for (int i = 0; i < n; ++i) {

        int current =
            path[i] - 1;

        int next =
            path[(i + 1) % n] - 1;


        // Solo la sfera usa la distanza geodetica
        if (geometry == Geometry::Sphere) {

            length +=
                GeodesicDistance(
                    cities[current],
                    cities[next]
                );
        }

        else {

            length +=
                EuclideanDistance(
                    cities[current],
                    cities[next]
                );
        }
    }

    return length;
}


// ============================================================
// Valutazione della popolazione
// ============================================================

void EvaluatePopulation(
    vector<Individual>& population,
    const vector<City>& cities,
    Geometry geometry
) {

    for (Individual& individual : population) {

        individual.length =
            PathLength(
                individual.path,
                cities,
                geometry
            );
    }


    // Ordinamento crescente di L
    sort(
        population.begin(),
        population.end(),
        [](const Individual& a,
           const Individual& b) {

            return a.length < b.length;
        }
    );
}


// ============================================================
// Media sulla migliore metà della popolazione
// ============================================================

double BestHalfMean(
    const vector<Individual>& population
) {

    int half =
        static_cast<int>(
            population.size()
        ) / 2;


    double sum = 0.0;


    for (int i = 0; i < half; ++i) {

        sum += population[i].length;
    }


    return
        sum / static_cast<double>(half);
}


// ============================================================
// Selezione
// ============================================================

int SelectIndividual(
    Random& rnd,
    int population_size,
    double selection_power
) {

    double r = rnd.Rannyu();


    int index =
        static_cast<int>(
            population_size *
            pow(r, selection_power)
        );


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

    int n =
        static_cast<int>(
            individual.path.size()
        );


    int i =
        RandomInt(rnd, 1, n - 1);

    int j =
        RandomInt(rnd, 1, n - 1);


    while (j == i) {

        j =
            RandomInt(rnd, 1, n - 1);
    }


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

    int n =
        static_cast<int>(
            individual.path.size()
        );


    int movable = n - 1;


    int m =
        RandomInt(
            rnd,
            1,
            movable - 1
        );


    int shift =
        RandomInt(
            rnd,
            1,
            movable - m
        );


    int total =
        m + shift;


    int start =
        RandomInt(
            rnd,
            1,
            n - total
        );


    rotate(
        individual.path.begin() + start,
        individual.path.begin() + start + m,
        individual.path.begin() + start + total
    );
}


// ============================================================
// MUTAZIONE 3:
// scambio di due blocchi
// ============================================================

void MutationBlockSwap(
    Individual& individual,
    Random& rnd
) {

    int n =
        static_cast<int>(
            individual.path.size()
        );


    int movable = n - 1;


    int m =
        RandomInt(
            rnd,
            1,
            movable / 2
        );


    int first_start =
        RandomInt(
            rnd,
            1,
            n - m
        );


    vector<int> possible_starts;


    for (int second_start = 1;
         second_start <= n - m;
         ++second_start) {

        bool non_overlapping =
            (second_start + m - 1 < first_start) ||
            (first_start + m - 1 < second_start);


        if (non_overlapping) {

            possible_starts.push_back(
                second_start
            );
        }
    }


    if (possible_starts.empty())
        return;


    int chosen =
        RandomInt(
            rnd,
            0,
            static_cast<int>(
                possible_starts.size()
            ) - 1
        );


    int second_start =
        possible_starts[chosen];


    for (int k = 0; k < m; ++k) {

        swap(
            individual.path[first_start + k],
            individual.path[second_start + k]
        );
    }
}


// ============================================================
// MUTAZIONE 4:
// inversione
// ============================================================

void MutationInversion(
    Individual& individual,
    Random& rnd
) {

    int n =
        static_cast<int>(
            individual.path.size()
        );


    int movable = n - 1;


    int m =
        RandomInt(
            rnd,
            2,
            movable
        );


    int start =
        RandomInt(
            rnd,
            1,
            n - m
        );


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

    if (rnd.Rannyu() < parameters.p_pair_swap) {

        MutationPairSwap(
            individual,
            rnd
        );
    }


    if (rnd.Rannyu() < parameters.p_shift) {

        MutationShift(
            individual,
            rnd
        );
    }


    if (rnd.Rannyu() < parameters.p_block_swap) {

        MutationBlockSwap(
            individual,
            rnd
        );
    }


    if (rnd.Rannyu() < parameters.p_inversion) {

        MutationInversion(
            individual,
            rnd
        );
    }
}


// ============================================================
// Crossover
// ============================================================

pair<Individual, Individual> Crossover(
    const Individual& parent1,
    const Individual& parent2,
    Random& rnd
) {

    int n =
        static_cast<int>(
            parent1.path.size()
        );


    int cut =
        RandomInt(
            rnd,
            1,
            n - 1
        );


    Individual child1;
    Individual child2;


    child1.path.assign(n, -1);
    child2.path.assign(n, -1);


    // ========================================================
    // Primo figlio
    // ========================================================

    vector<bool> used1(
        n + 1,
        false
    );


    for (int i = 0; i < cut; ++i) {

        child1.path[i] =
            parent1.path[i];

        used1[parent1.path[i]] =
            true;
    }


    int position = cut;


    for (int city : parent2.path) {

        if (!used1[city]) {

            child1.path[position] =
                city;

            used1[city] = true;

            position++;
        }
    }


    // ========================================================
    // Secondo figlio
    // ========================================================

    vector<bool> used2(
        n + 1,
        false
    );


    for (int i = 0; i < cut; ++i) {

        child2.path[i] =
            parent2.path[i];

        used2[parent2.path[i]] =
            true;
    }


    position = cut;


    for (int city : parent1.path) {

        if (!used2[city]) {

            child2.path[position] =
                city;

            used2[city] = true;

            position++;
        }
    }


    return {
        child1,
        child2
    };
}


// ============================================================
// Popolazione iniziale
// ============================================================

vector<Individual> CreateInitialPopulation(
    Random& rnd,
    const GAParameters& parameters
) {

    vector<Individual> population;


    population.reserve(
        parameters.population_size
    );


    for (int i = 0;
         i < parameters.population_size;
         ++i) {

        Individual individual;


        individual.path =
            RandomPath(
                rnd,
                parameters.n_cities
            );


        if (
            !CheckPath(
                individual.path,
                parameters.n_cities
            )
        ) {

            cerr
                << "ERRORE: cromosoma iniziale non valido."
                << endl;

            exit(1);
        }


        population.push_back(
            individual
        );
    }


    return population;
}


// ============================================================
// Generazione successiva
// ============================================================

vector<Individual> NewGeneration(
    const vector<Individual>& population,
    Random& rnd,
    const GAParameters& parameters
) {

    vector<Individual> new_population;


    new_population.reserve(
        parameters.population_size
    );


    // ========================================================
    // Elitismo
    // ========================================================

    for (int i = 0;
         i < parameters.n_elite &&
         i < parameters.population_size;
         ++i) {

        new_population.push_back(
            population[i]
        );
    }


    // ========================================================
    // Produzione dei figli
    // ========================================================

    while (
        static_cast<int>(
            new_population.size()
        )
        < parameters.population_size
    ) {

        // Selezione
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


        const Individual& parent1 =
            population[index1];

        const Individual& parent2 =
            population[index2];


        Individual child1;
        Individual child2;


        // Crossover
        if (
            rnd.Rannyu()
            < parameters.p_crossover
        ) {

            auto children =
                Crossover(
                    parent1,
                    parent2,
                    rnd
                );


            child1 =
                children.first;

            child2 =
                children.second;
        }

        else {

            child1 = parent1;
            child2 = parent2;
        }


        // Mutazioni
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


        // Controllo validità
        if (
            !CheckPath(
                child1.path,
                parameters.n_cities
            )
        ) {

            cerr
                << "ERRORE: primo figlio non valido."
                << endl;

            exit(1);
        }


        if (
            !CheckPath(
                child2.path,
                parameters.n_cities
            )
        ) {

            cerr
                << "ERRORE: secondo figlio non valido."
                << endl;

            exit(1);
        }


        new_population.push_back(
            child1
        );


        if (
            static_cast<int>(
                new_population.size()
            )
            < parameters.population_size
        ) {

            new_population.push_back(
                child2
            );
        }
    }


    return new_population;
}


// ============================================================
// Controllo delle generazioni snapshot
// ============================================================

bool IsSnapshotGeneration(
    int generation,
    const vector<int>& snapshot_generations
) {

    return
        find(
            snapshot_generations.begin(),
            snapshot_generations.end(),
            generation
        )
        != snapshot_generations.end();
}


// ============================================================
// Salvataggio delle città
// ============================================================

void SaveCities(
    const vector<City>& cities,
    const string& filename,
    Geometry geometry
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
        << scientific
        << setprecision(12);


    // Sfera e cubo sono configurazioni tridimensionali
    if (
        geometry == Geometry::Sphere ||
        geometry == Geometry::Cube
    ) {

        output
            << "# Coordinate delle citta' in 3D\n";

        output << "#\n";
        output << "# city x y z\n";


        for (size_t i = 0;
             i < cities.size();
             ++i) {

            output
                << i + 1 << " "
                << cities[i].x << " "
                << cities[i].y << " "
                << cities[i].z
                << "\n";
        }
    }

    else {

        output
            << "# Coordinate delle citta'\n";

        output << "#\n";
        output << "# city x y\n";


        for (size_t i = 0;
             i < cities.size();
             ++i) {

            output
                << i + 1 << " "
                << cities[i].x << " "
                << cities[i].y
                << "\n";
        }
    }


    output.close();
}


// ============================================================
// Salvataggio del miglior percorso
// ============================================================

void SaveBestPath(
    const Individual& best,
    const vector<City>& cities,
    const string& filename,
    Geometry geometry
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


    bool is_3D =
        geometry == Geometry::Sphere ||
        geometry == Geometry::Cube;


    if (is_3D) {

        output
            << "# order city x y z\n";
    }

    else {

        output
            << "# order city x y\n";
    }


    output
        << scientific
        << setprecision(12);


    for (size_t i = 0;
         i < best.path.size();
         ++i) {

        int city =
            best.path[i];


        output
            << i + 1 << " "
            << city << " "
            << cities[city - 1].x << " "
            << cities[city - 1].y;


        if (is_3D) {

            output
                << " "
                << cities[city - 1].z;
        }


        output << "\n";
    }


    // Chiusura del percorso
    int first_city =
        best.path[0];


    output
        << best.path.size() + 1
        << " "
        << first_city
        << " "
        << cities[first_city - 1].x
        << " "
        << cities[first_city - 1].y;


    if (is_3D) {

        output
            << " "
            << cities[first_city - 1].z;
    }


    output << "\n";


    output.close();
}


// ============================================================
// Esecuzione del Genetic Algorithm
// ============================================================

void RunGeneticAlgorithm(
    Random& rnd,
    const vector<City>& cities,
    const GAParameters& parameters,
    const string& name,
    Geometry geometry
) {

    cout << "\n============================================\n";
    cout << "Simulazione: " << name << endl;
    cout << "============================================\n";


    // Popolazione iniziale
    vector<Individual> population =
        CreateInitialPopulation(
            rnd,
            parameters
        );


    EvaluatePopulation(
        population,
        cities,
        geometry
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


    if (
        IsSnapshotGeneration(
            0,
            parameters.snapshot_generations
        )
    ) {

        SaveBestPath(
            population[0],
            cities,
            name + "_best_path_gen_0.txt",
            geometry
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

        population =
            NewGeneration(
                population,
                rnd,
                parameters
            );


        EvaluatePopulation(
            population,
            cities,
            geometry
        );


        double best_length =
            population[0].length;


        double best_half_mean =
            BestHalfMean(
                population
            );


        convergence
            << generation << " "
            << best_length << " "
            << best_half_mean
            << "\n";


        // Snapshot del miglior individuo
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
                filename,
                geometry
            );
        }


        // Output a terminale ogni 100 generazioni
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
    // Miglior percorso finale
    // ========================================================

    SaveBestPath(
        population[0],
        cities,
        name + "_best_path.txt",
        geometry
    );


    cout
        << "\nL_best finale = "
        << population[0].length
        << endl;


    cout << "Percorso:\n";


    for (int city : population[0].path) {

        cout
            << city
            << " -> ";
    }


    cout
        << population[0].path[0]
        << endl;
}


// ============================================================
// MAIN
// ============================================================

int main() {

    // ========================================================
    // Generatore random
    // ========================================================

    Random rnd;

    InitializeRandom(rnd);


    // ========================================================
    // Parametri del Genetic Algorithm
    // ========================================================

    GAParameters parameters;


    // ========================================================
    // Generazione di tutte le configurazioni
    //
    // Le città vengono generate prima di avviare i GA,
    // in modo che la loro geometria non dipenda dal numero
    // di estrazioni effettuate durante le simulazioni precedenti.
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


    vector<City> sphere_cities =
        GenerateCitiesSphere(
            rnd,
            parameters.n_cities
        );


    vector<City> cube_cities =
        GenerateCitiesCube(
            rnd,
            parameters.n_cities
        );


    // ========================================================
    // Salvataggio delle città
    // ========================================================

    SaveCities(
        circle_cities,
        "circle_cities.txt",
        Geometry::Circle
    );


    SaveCities(
        square_cities,
        "square_cities.txt",
        Geometry::Square
    );


    SaveCities(
        sphere_cities,
        "sphere_cities.txt",
        Geometry::Sphere
    );


    SaveCities(
        cube_cities,
        "cube_cities.txt",
        Geometry::Cube
    );


    // ========================================================
    // Genetic Algorithm - circonferenza
    // ========================================================

    RunGeneticAlgorithm(
        rnd,
        circle_cities,
        parameters,
        "circle",
        Geometry::Circle
    );


    // ========================================================
    // Genetic Algorithm - quadrato
    // ========================================================

    RunGeneticAlgorithm(
        rnd,
        square_cities,
        parameters,
        "square",
        Geometry::Square
    );


    // ========================================================
    // Genetic Algorithm - sfera
    // ========================================================

    RunGeneticAlgorithm(
        rnd,
        sphere_cities,
        parameters,
        "sphere",
        Geometry::Sphere
    );


    // ========================================================
    // Genetic Algorithm - cubo
    // ========================================================

    RunGeneticAlgorithm(
        rnd,
        cube_cities,
        parameters,
        "cube",
        Geometry::Cube
    );


    // ========================================================
    // Salvataggio del seed finale
    // ========================================================

    rnd.SaveSeed();


    cout << "\n============================================\n";
    cout << "Tutte le simulazioni sono state completate." << endl;
    cout << "============================================\n";


    return 0;
}