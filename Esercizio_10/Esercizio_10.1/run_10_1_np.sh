#!/usr/bin/env bash
# =============================================================
# Laboratorio di Simulazione Numerica - Esercizio 10.1
# Esegue le simulazioni su circonferenza e quadrato con p rank.
# Salva i dati separatamente dai vecchi test con 4 processi.
# =============================================================
set -euo pipefail

NPROC="${NPROC:-10}"
GENERATIONS="${GENERATIONS:-1000}"
POPULATION="${POPULATION:-500}"
INTERVAL="${INTERVAL:-50}"
TRIAL="${TRIAL:-0}"
FORCE="${FORCE:-0}"

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="$ROOT_DIR/results_10_1_np${NPROC}"

for number in "$NPROC" "$GENERATIONS" "$POPULATION" "$INTERVAL"; do
    if ! [[ "$number" =~ ^[0-9]+$ ]] || (( number < 1 )); then
        echo "ERRORE: numero di processi, generazioni, popolazione e intervallo devono essere positivi." >&2
        exit 1
    fi
done
if ! [[ "$TRIAL" =~ ^[0-9]+$ ]] || (( NPROC > 11 || POPULATION < 2 )); then
    echo "ERRORE: trial non valido oppure numero processi > 11." >&2
    exit 1
fi
for filename in main.exe Primes seed.in; do
    [[ -f "$ROOT_DIR/$filename" ]] || { echo "ERRORE: manca $filename" >&2; exit 1; }
done
command -v mpiexec >/dev/null 2>&1 || { echo "ERRORE: mpiexec non trovato" >&2; exit 1; }

mkdir -p "$RESULTS_DIR"
SETTINGS_FILE="$RESULTS_DIR/experiment_settings.txt"
SETTINGS="nprocesses=$NPROC
generations=$GENERATIONS
population_per_rank=$POPULATION
interval=$INTERVAL
trial=$TRIAL
geometry=both"
if [[ -f "$SETTINGS_FILE" && "$(cat "$SETTINGS_FILE")" != "$SETTINGS" ]]; then
    echo "ERRORE: la cartella $RESULTS_DIR contiene esperimenti con parametri diversi." >&2
    echo "Rinomina la cartella prima di modificarli." >&2
    exit 1
fi
printf '%s\n' "$SETTINGS" > "$SETTINGS_FILE"

for mode in independent migration; do
    CASE_DIR="$RESULTS_DIR/$mode"
    mkdir -p "$CASE_DIR"
    cp "$ROOT_DIR/main.exe" "$ROOT_DIR/Primes" "$ROOT_DIR/seed.in" "$CASE_DIR/"
    LOGFILE="$CASE_DIR/console.log"
    if [[ "$FORCE" != "1" && -f "$LOGFILE" ]] && \
       [[ -f "$CASE_DIR/circle_${mode}_trial_${TRIAL}_global_convergence.txt" ]] && \
       [[ -f "$CASE_DIR/square_${mode}_trial_${TRIAL}_global_convergence.txt" ]] && \
       [[ "$(grep -c 'RISULTATO FINALE' "$LOGFILE")" -eq 2 ]]; then
        echo "[SKIP] 10.1 | $mode (gia' completato)"
        continue
    fi
    echo "[RUN ] 10.1 | $mode | $NPROC processi | $GENERATIONS generazioni"
    if ! (
        cd "$CASE_DIR"
        mpiexec -np "$NPROC" ./main.exe \
            --geometry both --mode "$mode" --interval "$INTERVAL" \
            --generations "$GENERATIONS" --population "$POPULATION" \
            --trial "$TRIAL" > console.log 2>&1
    ); then
        echo "ERRORE: controlla $LOGFILE" >&2
        exit 1
    fi
    grep -E 'RISULTATO FINALE|L_best globale:|Tempo:' "$LOGFILE" || true
    [[ "$(grep -c 'RISULTATO FINALE' "$LOGFILE")" -eq 2 ]] || {
        echo "ERRORE: output incompleto in $LOGFILE" >&2
        exit 1
    }
done

echo
echo "Esercizio 10.1 completato: $RESULTS_DIR"
