#!/usr/bin/env python3

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt


def read_energy_log(csv_path: Path) -> tuple[list[float], list[float]]:
    times: list[float] = []
    energies: list[float] = []

    with csv_path.open(newline="") as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            times.append(float(row["time"]))
            energies.append(float(row["energy"]))

    return times, energies


def main() -> int:
    if len(sys.argv) != 3:
        print(
            f"Usage: {sys.argv[0]} <energy_log.csv> <output.png>",
            file=sys.stderr,
        )
        return 2

    csv_path = Path(sys.argv[1])
    output_path = Path(sys.argv[2])
    output_path.parent.mkdir(parents=True, exist_ok=True)

    times, energies = read_energy_log(csv_path)

    # plt.style.use("dark_background")
    fig, ax = plt.subplots(figsize=(10, 6), dpi=180)
    # fig.patch.set_facecolor("#111111")
    # ax.set_facecolor("#111111")
    ax.plot(times, energies,  linewidth=2.0)
    ax.set_xlabel("time")
    ax.set_ylabel("energy")
    ax.set_title("Energy over Time")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())