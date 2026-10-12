#!/usr/bin/env bash
# ============================================================
# Laboratorio di Simulazione Numerica - Esercizio 10.2
# Benchmark MPI: GA indipendenti vs GA con migrazioni
# Da avviare dal terminale MSYS2 MINGW64 con:
#     bash run_trials.sh
# ============================================================

set -euo pipefail

# Parametri del confronto. Modificabili all'avvio, ad esempio:
# TRIALS=3 GENERATIONS=50 bash run_trials.sh
TRIALS="${TRIALS:-20}"
NPROC="${NPROC:-4}"
GENERATIONS="${GENERATIONS:-1000}"
POPULATION="${POPULATION:-500}"
FORCE="${FORCE:-0}"  # FORCE=1 ripete anche le prove completate

# Intervalli di migrazione da confrontare.
INTERVALS=(10 20 50 100)

# La cartella di lavoro e' quella in cui si trova lo script.
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="$ROOT_DIR/results_10_2"

# Controllo dei parametri.
for number in "$TRIALS" "$NPROC" "$GENERATIONS" "$POPULATION"; do
    if ! [[ "$number" =~ ^[0-9]+$ ]] || (( number < 1 )); then
        echo "ERRORE: tutti i parametri numerici devono essere interi positivi." >&2
        exit 1
    fi
done

if (( NPROC > 11 )); then
    echo "ERRORE: il programma C++ supporta al massimo 11 processi." >&2
    exit 1
fi

if (( POPULATION < 2 )); then
    echo "ERRORE: servono almeno 2 individui per popolazione." >&2
    exit 1
fi

if ! command -v mpiexec >/dev/null 2>&1; then
    echo "ERRORE: mpiexec non trovato. Controlla il PATH di MS MPI." >&2
    exit 1
fi

for filename in main.exe Primes seed.in cap_prov_ita.dat; do
    if [[ ! -f "$ROOT_DIR/$filename" ]]; then
        echo "ERRORE: manca $filename nella cartella $ROOT_DIR" >&2
        exit 1
    fi
done

mkdir -p "$RESULTS_DIR"

# Non mischiare output di esperimenti con parametri differenti.
SETTINGS_FILE="$RESULTS_DIR/experiment_settings.txt"
if [[ -f "$SETTINGS_FILE" && "$FORCE" != "1" ]]; then
    for entry in "nprocesses=$NPROC" "ntrials=$TRIALS" \
                 "generations=$GENERATIONS" "population_per_rank=$POPULATION"; do
        if ! grep -Fxq "$entry" "$SETTINGS_FILE"; then
            echo "ERRORE: esistono risultati con parametri diversi in $RESULTS_DIR" >&2
            echo "Per ricominciare il benchmark con nuovi parametri usa FORCE=1." >&2
            exit 1
        fi
    done
fi

# Salva i parametri comuni in formato facilmente consultabile.
{
    echo "nprocesses=$NPROC"
    echo "ntrials=$TRIALS"
    echo "generations=$GENERATIONS"
    echo "population_per_rank=$POPULATION"
    echo "intervals=${INTERVALS[*]}"
    echo "geometry=italy"
    echo "note=Confronti appaiati per trial; la stessa configurazione iniziale e' usata per tutte le modalita'"
} > "$SETTINGS_FILE"

# Ogni configurazione ha la propria directory, perche' il programma
# originale non include --interval nei nomi dei file di output.
run_case() {
    local mode="$1"
    local interval="$2"
    local label="$3"
    local case_dir="$RESULTS_DIR/$label"
    mkdir -p "$case_dir"

    # I processi Windows lanciati da mpiexec lavorano dentro case_dir.
    # Devono quindi trovare anche i file di input nella stessa cartella.
    cp "$ROOT_DIR/main.exe" "$case_dir/main.exe"
    cp "$ROOT_DIR/Primes" "$ROOT_DIR/seed.in" \
       "$ROOT_DIR/cap_prov_ita.dat" "$case_dir/"

    for ((trial=0; trial<TRIALS; trial++)); do
        local prefix="italy_${mode}_trial_${trial}"
        local logfile="$case_dir/${prefix}_console.log"
        local convfile="$case_dir/${prefix}_global_convergence.txt"

        if [[ "$FORCE" != "1" && -f "$logfile" && -f "$convfile" ]] && \
           grep -q "RISULTATO FINALE" "$logfile"; then
            echo "[SKIP] $label | trial $trial gia' completato"
            continue
        fi

        echo "[RUN ] $label | trial $trial / $((TRIALS-1))"
        # Subshell: il cd non modifica la cartella di chi esegue lo script.
        # Sono usati gli stessi trial per tutti gli intervalli.
        if ! (
            cd "$case_dir"
            mpiexec -np "$NPROC" ./main.exe \
                --geometry italy \
                --mode "$mode" \
                --interval "$interval" \
                --generations "$GENERATIONS" \
                --population "$POPULATION" \
                --trial "$trial" \
                > "${prefix}_console.log" 2>&1
        ); then
            echo "ERRORE durante $label, trial $trial." >&2
            echo "Controlla: $logfile" >&2
            exit 1
        fi

        if [[ ! -f "$convfile" ]] || ! grep -q "RISULTATO FINALE" "$logfile"; then
            echo "ERRORE: output mancante o incompleto: $logfile" >&2
            exit 1
        fi
    done
}

# 20 trial senza migrazioni: riferimento comune a tutti i confronti.
run_case "independent" "50" "independent"

# 20 trial per ciascuna frequenza di migrazione.
for interval in "${INTERVALS[@]}"; do
    run_case "migration" "$interval" "migration_${interval}"
done

echo
echo "Tutti gli esperimenti sono completati."
echo "File salvati in: $RESULTS_DIR"
echo "Per creare le tabelle di confronto:"
echo "    python analyze_trials.py"
