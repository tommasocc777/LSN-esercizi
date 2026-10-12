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
#include <mpi.h>
#include <limits>

using namespace std;

// ============================================================

// Strutture dati

// ============================================================

// Coordinate x,y (cerchio/quadrato) oppure longitudine,latitudine (Italia)

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

void InitializeRandom(Random& rnd, int stream, int trial = 0) {

    // Ogni stream utilizza una diversa coppia di numeri di Primes.
    // Lo stesso seed.in viene mantenuto per la riproducibilita'.
    int seed[4] = {0, 0, 0, 1};
    int p1 = 0, p2 = 0;

    ifstream primes("Primes");
    if (!primes.is_open()) {
        cerr << "ERRORE: impossibile aprire Primes." << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Saltiamo le coppie precedenti e leggiamo la coppia scelta.
    for (int i = 0; i <= stream; ++i) {
        if (!(primes >> p1 >> p2)) {
            cerr << "ERRORE: Primes non contiene abbastanza coppie."
                 << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    ifstream input("seed.in");
    if (!input.is_open()) {
        cerr << "ERRORE: impossibile aprire seed.in." << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    bool found = false;
    string property;
    while (input >> property) {
        if (property == "RANDOMSEED") {
            if (!(input >> seed[0] >> seed[1] >> seed[2] >> seed[3])) {
                cerr << "ERRORE: valori RANDOMSEED non validi." << endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            found = true;
            break;
        }
    }

    if (!found) {
        cerr << "ERRORE: RANDOMSEED non trovato in seed.in." << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Trial differenti, riproducibili: a parita' di trial e rank, le
    // modalita' independent e migration partono dalla stessa popolazione.
    // Le sequenze non sono matematicamente garantite indipendenti.
    seed[0] = (seed[0] + 37 * (trial % 4096)) % 4096;
    seed[1] = (seed[1] + 101 * (trial % 4096)) % 4096;
    seed[2] = (seed[2] + 229 * (trial % 4096)) % 4096;
    seed[3] = (seed[3] + 487 * (trial % 4096)) % 4096;

    rnd.SetRandom(seed, p1, p2);
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

// ============================================================
// Lettura delle 110 coordinate geografiche
// ============================================================

vector<City> ReadItalianCities(const string& filename) {

    ifstream input(filename);
    if (!input.is_open()) {
        cerr << "ERRORE: impossibile aprire " << filename << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    vector<City> cities;
    City city;

    // Nel file: prima LONGITUDINE (x), poi LATITUDINE (y).
    while (input >> city.x >> city.y) {
        if (!isfinite(city.x) || !isfinite(city.y) ||
            city.x < -180.0 || city.x > 180.0 ||
            city.y < -90.0 || city.y > 90.0) {
            cerr << "ERRORE: coordinate geografiche non valide." << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        cities.push_back(city);
    }

    // Ogni riga deve contenere esattamente due numeri validi.
    if (!input.eof() || cities.size() != 110) {
        cerr << "ERRORE: attese 110 coppie lon-lat valide in "
             << filename << "." << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    return cities;
}

// ============================================================
// Matrice delle distanze: coordinate cartesiane o geografiche
// ============================================================

using DistanceMatrix = vector<vector<double>>;

double GeographicDistanceKm(const City& a, const City& b) {

    // Haversine: distanza ortodromica su Terra sferica.
    const double R = 6371.0;
    const double DEG_TO_RAD = acos(-1.0) / 180.0;

    double phi1 = a.y * DEG_TO_RAD;
    double phi2 = b.y * DEG_TO_RAD;
    double dphi = (b.y - a.y) * DEG_TO_RAD;
    double dlambda = (b.x - a.x) * DEG_TO_RAD;

    double h = pow(sin(dphi / 2.0), 2)
             + cos(phi1) * cos(phi2) * pow(sin(dlambda / 2.0), 2);

    // Protezione da piccoli errori di arrotondamento.
    h = max(0.0, min(1.0, h));

    return 2.0 * R * asin(sqrt(h));
}

DistanceMatrix BuildDistanceMatrix(
    const vector<City>& cities,
    bool geographic
) {

    int n = static_cast<int>(cities.size());
    DistanceMatrix distances(n, vector<double>(n, 0.0));

    // Calcoliamo ciascuna distanza una sola volta. Le distanze
    // rimangono costanti durante tutte le generazioni del GA.
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double distance;

            if (geographic) {
                distance = GeographicDistanceKm(cities[i], cities[j]);
            } else {
                double dx = cities[i].x - cities[j].x;
                double dy = cities[i].y - cities[j].y;
                distance = sqrt(dx * dx + dy * dy);
            }

            distances[i][j] = distance;
            distances[j][i] = distance;
        }
    }

    return distances;
}

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
    const DistanceMatrix& distances
) {

    double length = 0.0;
    int n_cities = static_cast<int>(path.size());

    for (int i = 0; i < n_cities; ++i) {

        // La citta' 1 e' sempre la prima; il modulo chiude il tour.
        int current = path[i] - 1;
        int next = path[(i + 1) % n_cities] - 1;

        length += distances[current][next];
    }

    return length;
}

// ============================================================

// Valutazione e ordinamento della popolazione

// ============================================================

void EvaluatePopulation(

    vector<Individual>& population,

    const DistanceMatrix& distances

) {

    // Calcolo della funzione costo di ogni individuo

    for (Individual& individual : population) {

        individual.length =

            PathLength(individual.path, distances);

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

    output << "# Colonna 2: x oppure longitudine (italy)\n";

    output << "# Colonna 3: y oppure latitudine (italy)\n";

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

    output << "# Colonna 3: x oppure longitudine (italy)\n";

    output << "# Colonna 4: y oppure latitudine (italy)\n";

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
// Geometrie comuni a tutti i processi MPI
// ============================================================

vector<City> BroadcastCities(
    const vector<City>& root_cities,
    int n_cities,
    int rank
) {

    vector<double> coordinates(2 * n_cities, 0.0);

    // Solo il rank 0 conosce inizialmente le coordinate.
    if (rank == 0) {
        for (int i = 0; i < n_cities; ++i) {
            coordinates[2 * i] = root_cities[i].x;
            coordinates[2 * i + 1] = root_cities[i].y;
        }
    }

    // Tutti ricevono le stesse coordinate.
    MPI_Bcast(coordinates.data(), 2 * n_cities, MPI_DOUBLE,
              0, MPI_COMM_WORLD);

    vector<City> cities(n_cities);
    for (int i = 0; i < n_cities; ++i) {
        cities[i].x = coordinates[2 * i];
        cities[i].y = coordinates[2 * i + 1];
    }

    return cities;
}


// ============================================================
// Migrazione tra continenti
// ============================================================

void Migrate(
    vector<Individual>& population,
    const DistanceMatrix& distances,
    Random& migration_rnd,
    int rank,
    int size
) {

    vector<int> partner(size, -1);

    if (rank == 0) {
        // Permutazione casuale dei rank: abbinamenti differenti
        // a ogni migrazione, senza fissare la topologia dei continenti.
        vector<int> order(size);
        iota(order.begin(), order.end(), 0);

        for (int i = size - 1; i > 0; --i) {
            int j = RandomInt(migration_rnd, 0, i);
            swap(order[i], order[j]);
        }

        // Per un numero dispari un continente rimane escluso
        // soltanto da questo turno di migrazione.
        for (int i = 0; i + 1 < size; i += 2) {
            partner[order[i]] = order[i + 1];
            partner[order[i + 1]] = order[i];
        }
    }

    MPI_Bcast(partner.data(), size, MPI_INT, 0, MPI_COMM_WORLD);

    int destination = partner[rank];
    if (destination < 0) return;

    // La popolazione e' ordinata per lunghezza crescente.
    // Scambiamo l'intero cromosoma del miglior individuo.
    int n_cities = static_cast<int>(distances.size());
    vector<int> received_path(n_cities);

    MPI_Sendrecv(
        population[0].path.data(), n_cities, MPI_INT,
        destination, 100,
        received_path.data(), n_cities, MPI_INT,
        destination, 100,
        MPI_COMM_WORLD, MPI_STATUS_IGNORE
    );

    if (!CheckPath(received_path, n_cities)) {
        cerr << "ERRORE: percorso ricevuto non valido al rank "
             << rank << endl;
        MPI_Abort(MPI_COMM_WORLD, 2);
    }

    // Il migrante rimpiazza l'individuo peggiore.
    Individual migrant;
    migrant.path = move(received_path);
    migrant.length = PathLength(migrant.path, distances);

    population.back() = move(migrant);

    // Ricalcoliamo e riordiniamo la popolazione.
    EvaluatePopulation(population, distances);
}


// ============================================================
// Esecuzione di un Genetic Algorithm su ciascun continente
// ============================================================

void RunParallelGA(
    const vector<City>& cities,
    const GAParameters& parameters,
    const string& geometry_name,
    int geometry_index,
    bool use_migration,
    int migration_interval,
    int rank,
    int size,
    int trial
) {

    // Lo stream del GA dipende solo da geometria e rank,
    // NON dalla presenza delle migrazioni: confronto riproducibile.
    Random rnd;
    InitializeRandom(rnd, 1 + geometry_index * size + rank, trial);

    // La matrice evita di ricalcolare trigonometria e distanze
    // a ogni valutazione del percorso TSP.
    DistanceMatrix distances = BuildDistanceMatrix(
        cities, geometry_name == "italy");

    // Random separato dedicato esclusivamente agli abbinamenti.
    Random migration_rnd;
    if (rank == 0) {
        InitializeRandom(migration_rnd, 1 + 3 * size, trial);
    }

    string mode = use_migration ? "migration" : "independent";
    string prefix = geometry_name + "_" + mode
                  + "_trial_" + to_string(trial);
    string rank_prefix = prefix + "_rank_" + to_string(rank);

    ofstream convergence(rank_prefix + "_convergence.txt");
    if (!convergence.is_open()) {
        cerr << "ERRORE: file convergenza, rank " << rank << endl;
        MPI_Abort(MPI_COMM_WORLD, 3);
    }

    convergence << "# generation best_current best_ever best_half_mean\n";
    convergence << scientific << setprecision(12);

    // Un solo file globale scritto dal processo 0.
    ofstream global_convergence;
    if (rank == 0) {
        global_convergence.open(prefix + "_global_convergence.txt");
        if (!global_convergence.is_open()) {
            cerr << "ERRORE: file globale." << endl;
            MPI_Abort(MPI_COMM_WORLD, 3);
        }
        global_convergence << "# generation global_best_ever\n";
        global_convergence << scientific << setprecision(12);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    vector<Individual> population =
        CreateInitialPopulation(rnd, parameters);
    EvaluatePopulation(population, distances);

    Individual best_ever = population.front();

    // Stesse operazioni di riduzione per il controllo indipendente
    // e per quello con migrazione.
    for (int generation = 0;
         generation <= parameters.generations; ++generation) {

        if (generation > 0) {
            population = NewGeneration(population, rnd, parameters);
            EvaluatePopulation(population, distances);

            if (use_migration && size > 1 &&
                generation % migration_interval == 0) {
                Migrate(population, distances, migration_rnd, rank, size);
            }

            if (population.front().length < best_ever.length) {
                best_ever = population.front();
            }
        }

        // Per n_elite=0 il best corrente puo' peggiorare;
        // il miglior percorso storico invece non va perso.
        convergence << generation << " "
                    << population.front().length << " "
                    << best_ever.length << " "
                    << BestHalfMean(population) << "\n";

        if (IsSnapshotGeneration(generation,
                                 parameters.snapshot_generations)) {
            SaveBestPath(best_ever, cities,
                         rank_prefix + "_best_path_gen_" +
                         to_string(generation) + ".txt");
        }

        // Riduzione occasionale, non a ogni generazione,
        // per limitare l'overhead delle comunicazioni.
        if (generation % 10 == 0 ||
            generation == parameters.generations) {
            double local_best = best_ever.length;
            double global_best = 0.0;

            MPI_Reduce(&local_best, &global_best, 1, MPI_DOUBLE,
                       MPI_MIN, 0, MPI_COMM_WORLD);

            if (rank == 0) {
                global_convergence << generation << " "
                                   << global_best << "\n";
                if (generation % 100 == 0 ||
                    generation == parameters.generations) {
                    cout << geometry_name << " - " << mode
                         << " - generazione " << setw(4) << generation
                         << " - L_best = " << global_best << endl;
                }
            }
        }
    }

    convergence.close();
    if (rank == 0) global_convergence.close();

    SaveBestPath(best_ever, cities,
                 rank_prefix + "_best_path.txt");

    // Selezioniamo il processo con il migliore percorso storico.
    // In caso di parita' MPI_MINLOC sceglie il rank piu' piccolo.
    struct {
        double length;
        int rank;
    } candidate, winner;

    candidate.length = best_ever.length;
    candidate.rank = rank;

    MPI_Allreduce(&candidate, &winner, 1, MPI_DOUBLE_INT,
                  MPI_MINLOC, MPI_COMM_WORLD);

    vector<int> global_path(parameters.n_cities);
    if (rank == winner.rank) {
        global_path = best_ever.path;
    }

    MPI_Bcast(global_path.data(), parameters.n_cities,
              MPI_INT, winner.rank, MPI_COMM_WORLD);

    if (rank == 0) {
        if (!CheckPath(global_path, parameters.n_cities)) {
            cerr << "ERRORE: migliore percorso globale non valido."
                 << endl;
            MPI_Abort(MPI_COMM_WORLD, 4);
        }

        Individual global_best;
        global_best.path = move(global_path);
        global_best.length = PathLength(global_best.path, distances);
        SaveBestPath(global_best, cities,
                     prefix + "_global_best_path.txt");
    }

    // Il tempo di una simulazione e' quello del rank piu' lento.
    MPI_Barrier(MPI_COMM_WORLD);
    double local_time = MPI_Wtime() - start_time;
    double total_time = 0.0;
    MPI_Reduce(&local_time, &total_time, 1,
               MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        cout << "\nRISULTATO FINALE: " << prefix << "\n"
             << "  Processi: " << size << "\n"
             << "  L_best globale: " << winner.length << "\n"
             << "  Continent migliore: " << winner.rank << "\n"
             << "  Tempo: " << fixed << setprecision(3)
             << total_time << " secondi\n\n";
    }
}


// ============================================================
// MAIN - opzioni di esecuzione
// ============================================================

int main(int argc, char* argv[]) {

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    GAParameters parameters;
    bool use_migration = true;
    int migration_interval = 50;
    string geometry = "italy";
    string data_file = "cap_prov_ita.dat";
    int trial = 0;

    // Configurazione: tutti i rank ricevono gli stessi argomenti.
    // --mode independent|migration
    // --interval 50
    // --geometry italy|circle|square|both
    // --file cap_prov_ita.dat
    // --trial 0
    // --generations 1000
    // --population 500
    try {
        for (int i = 1; i < argc; ++i) {
            string option = argv[i];
            if (i + 1 >= argc) throw invalid_argument("Argomento mancante: " + option);
            string value = argv[++i];

            if (option == "--mode") {
                if (value == "migration") use_migration = true;
                else if (value == "independent") use_migration = false;
                else throw invalid_argument("--mode non valido");
            } else if (option == "--interval") {
                migration_interval = stoi(value);
            } else if (option == "--geometry") {
                geometry = value;
            } else if (option == "--file") {
                data_file = value;
            } else if (option == "--trial") {
                trial = stoi(value);
            } else if (option == "--generations") {
                parameters.generations = stoi(value);
            } else if (option == "--population") {
                parameters.population_size = stoi(value);
            } else {
                throw invalid_argument("Opzione sconosciuta: " + option);
            }
        }

        if (size < 1 || size > 11 || migration_interval <= 0 ||
            parameters.generations < 1 ||
            parameters.population_size < 2 || trial < 0 || trial > 10000 ||
            (geometry != "circle" && geometry != "square" &&
             geometry != "both" && geometry != "italy")) {
            throw invalid_argument("Parametri fuori intervallo");
        }
    } catch (const exception& err) {
        if (rank == 0) cerr << "ERRORE: " << err.what() << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Italia: file fisso di 110 coppie longitudine-latitudine.
    // Le geometrie sintetiche dell'Esercizio 9 rimangono disponibili.
    if (geometry == "italy") {

        parameters.n_cities = 110;
        vector<City> root_italy;

        if (rank == 0) {
            root_italy = ReadItalianCities(data_file);
        }

        vector<City> italy_cities =
            BroadcastCities(root_italy, parameters.n_cities, rank);

        if (rank == 0) {
            SaveCities(italy_cities, "italy_cities.txt");
        }

        RunParallelGA(italy_cities, parameters, "italy", 2,
                      use_migration, migration_interval, rank, size, trial);

    } else {

        // I punti sintetici sono sempre gli stessi in tutti i trial.
        Random geometry_rnd;
        vector<City> root_circle, root_square;

        if (rank == 0) {
            InitializeRandom(geometry_rnd, 0, 0);
            root_circle = GenerateCitiesCircle(geometry_rnd,
                                               parameters.n_cities);
            root_square = GenerateCitiesSquare(geometry_rnd,
                                               parameters.n_cities);
        }

        vector<City> circle_cities =
            BroadcastCities(root_circle, parameters.n_cities, rank);
        vector<City> square_cities =
            BroadcastCities(root_square, parameters.n_cities, rank);

        if (rank == 0) {
            SaveCities(circle_cities, "circle_cities.txt");
            SaveCities(square_cities, "square_cities.txt");
        }

        if (geometry == "circle" || geometry == "both") {
            RunParallelGA(circle_cities, parameters, "circle", 0,
                          use_migration, migration_interval, rank, size, trial);
        }

        if (geometry == "square" || geometry == "both") {
            RunParallelGA(square_cities, parameters, "square", 1,
                          use_migration, migration_interval, rank, size, trial);
        }
    }

    MPI_Finalize();
    return 0;
}
