#!/usr/bin/env python3
"""Assign explicit electrons to nearest carbon hosts under orthorhombic PBC.

This is the first analysis primitive needed for Jamil Fig. 3-style electron
mobility/path-density checks.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Atom:
    atom_id: int
    element: str
    x: float
    y: float
    z: float


def parse_geo(path: Path) -> tuple[list[Atom], tuple[float, float, float] | None]:
    atoms: list[Atom] = []
    box: tuple[float, float, float] | None = None
    with path.open("r", encoding="utf-8", errors="replace") as handle:
      for line in handle:
        if line.startswith("CRYSTX"):
            parts = line.split()
            if len(parts) >= 4:
                box = (float(parts[1]), float(parts[2]), float(parts[3]))
        elif line.startswith("HETATM"):
            parts = line.split()
            if len(parts) < 7:
                raise ValueError(f"Malformed HETATM line in {path}: {line.rstrip()}")
            atom_id = int(parts[1])
            element = parts[2]
            x, y, z = (float(parts[4]), float(parts[5]), float(parts[6]))
            atoms.append(Atom(atom_id, element, x, y, z))
    if not atoms:
        raise ValueError(f"No HETATM records found in {path}")
    return atoms, box


def minimum_image_delta(delta: float, box_length: float | None) -> float:
    if box_length is None or box_length <= 0.0:
        return delta
    return delta - box_length * round(delta / box_length)


def distance(a: Atom, b: Atom, box: tuple[float, float, float] | None) -> float:
    lengths = box if box is not None else (None, None, None)
    dx = minimum_image_delta(a.x - b.x, lengths[0])
    dy = minimum_image_delta(a.y - b.y, lengths[1])
    dz = minimum_image_delta(a.z - b.z, lengths[2])
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def assign_hosts(atoms: list[Atom], box: tuple[float, float, float] | None) -> list[dict[str, object]]:
    carbons = [atom for atom in atoms if atom.element == "C"]
    electrons = [atom for atom in atoms if atom.element == "El"]
    if not carbons:
        raise ValueError("No carbon atoms found")
    if not electrons:
        raise ValueError("No explicit electrons found")

    rows: list[dict[str, object]] = []
    for electron in electrons:
        host = min(carbons, key=lambda carbon: distance(electron, carbon, box))
        rows.append(
            {
                "electron_id": electron.atom_id,
                "host_carbon_id": host.atom_id,
                "distance_A": distance(electron, host, box),
                "electron_x": electron.x,
                "electron_y": electron.y,
                "electron_z": electron.z,
                "host_x": host.x,
                "host_y": host.y,
                "host_z": host.z,
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("geo", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    atoms, box = parse_geo(args.geo)
    rows = assign_hosts(atoms, box)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    unique_hosts = len({row["host_carbon_id"] for row in rows})
    print(f"electrons={len(rows)} unique_hosts={unique_hosts} output={args.out}")


if __name__ == "__main__":
    main()
