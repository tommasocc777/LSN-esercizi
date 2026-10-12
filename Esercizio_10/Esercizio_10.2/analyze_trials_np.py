#!/usr/bin/env python3
"""Esercizio 10.2 - analisi descrittiva dei trial MPI, senza librerie esterne.

Legge i file creati da run_trials.sh e salva:
  results_10_2/trial_results.csv
  results_10_2/paired_comparison.csv
  results_10_2/summary.csv

Eseguire nella cartella del progetto: python analyze_trials.py
"""

from __future__ import annotations

import argparse
import csv
import re
import statistics as stats
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results_10_2_np10"  # Impostato da --nproc nel main.


def read_settings() -> tuple[int, tuple[int, ...]]:
    file = RESULTS / "experiment_settings.txt"
    if not file.is_file():
        raise FileNotFoundError("Non trovo experiment_settings.txt: esegui prima run_trials.sh")
    settings = {}
    for line in file.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            settings[key] = value
    return int(settings["ntrials"]), tuple(map(int, settings["intervals"].split()))


def read_final_length(path: Path) -> float:
    """Legge l'ultimo valore non commentato della convergenza globale."""
    if not path.is_file():
        raise FileNotFoundError(f"Manca il file di convergenza: {path}")
    rows = [line.split() for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.lstrip().startswith("#")]
    if not rows:
        raise ValueError(f"File di convergenza vuoto: {path}")
    return float(rows[-1][1])


def read_console(path: Path) -> tuple[float, int]:
    """Estrae il tempo MPI in secondi e il rank migliore."""
    if not path.is_file():
        raise FileNotFoundError(f"Manca il log: {path}")
    content = path.read_text(encoding="utf-8", errors="replace")
    times = re.findall(r"Tempo:\s*([0-9]+(?:\.[0-9]+)?)\s+secondi", content)
    ranks = re.findall(r"Continent(?:e)? migliore:\s*(\d+)", content)
    if not times or not ranks or "RISULTATO FINALE" not in content:
        raise ValueError(f"Log incompleto: {path}")
    return float(times[-1]), int(ranks[-1])


def load_row(mode: str, interval: int | None, trial: int) -> dict:
    dirname = "independent" if mode == "independent" else f"migration_{interval}"
    folder = RESULTS / dirname
    prefix = f"italy_{mode}_trial_{trial}"
    length = read_final_length(folder / f"{prefix}_global_convergence.txt")
    elapsed, best_rank = read_console(folder / f"{prefix}_console.log")
    return {
        "mode": mode,
        "interval": "" if interval is None else interval,
        "trial": trial,
        "best_km": length,
        "time_s": elapsed,
        "best_rank": best_rank,
    }


def write_csv(filename: str, rows: list[dict], columns: list[str]) -> None:
    with (RESULTS / filename).open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def sample_std(values: list[float]) -> float:
    return stats.stdev(values) if len(values) > 1 else 0.0


def main() -> None:
    global RESULTS
    parser = argparse.ArgumentParser(description="Analisi benchmark MPI con p processi")
    parser.add_argument("--nproc", type=int, default=10, choices=range(1, 12),
                        metavar="1..11", help="Numero di processi MPI (default: 10)")
    args = parser.parse_args()
    RESULTS = ROOT / f"results_10_2_np{args.nproc}"
    settings_text = (RESULTS / "experiment_settings.txt").read_text(encoding="utf-8")
    if f"nprocesses={args.nproc}\n" not in (settings_text + "\n"):
        raise ValueError("Numero di processi non coerente con experiment_settings.txt")
    n_trials, intervals = read_settings()
    raw_rows = []
    control = {}
    migrated = {}

    for trial in range(n_trials):
        result = load_row("independent", None, trial)
        raw_rows.append(result)
        control[trial] = result

    for interval in intervals:
        for trial in range(n_trials):
            result = load_row("migration", interval, trial)
            raw_rows.append(result)
            migrated[(interval, trial)] = result

    write_csv("trial_results.csv", raw_rows,
              ["mode", "interval", "trial", "best_km", "time_s", "best_rank"])

    paired = []
    for interval in intervals:
        for trial in range(n_trials):
            baseline = control[trial]
            experiment = migrated[(interval, trial)]
            delta = baseline["best_km"] - experiment["best_km"]
            paired.append({
                "trial": trial,
                "interval": interval,
                "independent_km": baseline["best_km"],
                "migration_km": experiment["best_km"],
                "improvement_km": delta,
                "improvement_pct": 100.0 * delta / baseline["best_km"],
                "independent_time_s": baseline["time_s"],
                "migration_time_s": experiment["time_s"],
            })
    write_csv("paired_comparison.csv", paired,
              ["trial", "interval", "independent_km", "migration_km",
               "improvement_km", "improvement_pct", "independent_time_s",
               "migration_time_s"])

    summary = []
    for label, rows, interval in [
        ("independent", list(control.values()), None),
        *[("migration", [migrated[(p, t)] for t in range(n_trials)], p)
          for p in intervals],
    ]:
        lengths = [r["best_km"] for r in rows]
        times = [r["time_s"] for r in rows]
        differences = ([control[r["trial"]]["best_km"] - r["best_km"] for r in rows]
                       if interval is not None else [])
        wins = (sum(d > 1e-9 for d in differences) / n_trials * 100
                if interval is not None else "")
        summary.append({
            "mode": label,
            "interval": "" if interval is None else interval,
            "n_trials": n_trials,
            "mean_best_km": stats.mean(lengths),
            "std_best_km": sample_std(lengths),
            "median_best_km": stats.median(lengths),
            "min_best_km": min(lengths),
            "mean_time_s": stats.mean(times),
            "std_time_s": sample_std(times),
            "win_rate_pct": wins,
            "mean_improvement_km": stats.mean(differences) if differences else "",
        })

    write_csv("summary.csv", summary,
              ["mode", "interval", "n_trials", "mean_best_km",
               "std_best_km", "median_best_km", "min_best_km",
               "mean_time_s", "std_time_s", "win_rate_pct",
               "mean_improvement_km"])

    print("\nESERCIZIO 10.2 — Confronto su", n_trials, "trial per configurazione\n")
    print(f"{'Modalita':<18} {'L_media (km)':>15} {'Dev.std (km)':>15} "
          f"{'t_medio (s)':>14} {'Vittorie (%)':>14}")
    print("-" * 81)
    for row in summary:
        label = ("independent" if row["interval"] == ""
                 else f"migration {row['interval']}")
        wr = "-" if row["win_rate_pct"] == "" else f"{row['win_rate_pct']:.1f}"
        print(f"{label:<18} {row['mean_best_km']:>15.2f} "
              f"{row['std_best_km']:>15.2f} "
              f"{row['mean_time_s']:>14.3f} {wr:>14}")
    print("\nFile generati: trial_results.csv, paired_comparison.csv, summary.csv")
    print("Nota: la media piu' bassa descrive questi trial, non prova da sola che sia un ottimo universale.")


if __name__ == "__main__":
    main()
