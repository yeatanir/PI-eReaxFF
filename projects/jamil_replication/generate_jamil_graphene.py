#!/usr/bin/env python3
"""Generate Jamil-style paired-electron graphene geometries."""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path


A_CC = 1.4210
Z_C = 10.0
Z_EL = 10.0
PAIR_RADIUS = 0.35
LZ = 20.0
UNIT = [
    (0.0, 0.0),
    (0.5, 1.0 / 6.0),
    (0.5, 3.0 / 6.0),
    (0.0, 4.0 / 6.0),
]


@dataclass(frozen=True)
class Geometry:
    carbon: list[tuple[float, float, float]]
    electrons: list[tuple[float, float, float]]
    lx: float
    ly: float
    lz: float


def build_c_positions(nx: int, ny: int) -> tuple[list[tuple[float, float, float]], float, float]:
    coords = []
    ax = math.sqrt(3.0) * A_CC
    ay = 3.0 * A_CC
    lx = nx * ax
    ly = ny * ay
    for ix in range(nx):
        for iy in range(ny):
            for fx, fy in UNIT:
                coords.append(((ix + fx) * ax % lx, (iy + fy) * ay % ly, Z_C))
    return coords, lx, ly


def paired_offset(idx: int) -> tuple[float, float, float]:
    theta = 2.0 * math.pi * ((idx - 1) % 12) / 12.0
    return PAIR_RADIUS * math.cos(theta), PAIR_RADIUS * math.sin(theta), 0.0


def build_geometry(nx: int, ny: int, excess: bool) -> Geometry:
    carbon, lx, ly = build_c_positions(nx, ny)
    electrons = []
    for idx, (x, y, z) in enumerate(carbon, start=1):
        dx, dy, dz = paired_offset(idx)
        electrons.append(((x + dx) % lx, (y + dy) % ly, z + dz))
    if excess:
        # Put the excess electron away from the origin but still in-plane.
        electrons.append((0.70 * lx, 0.40 * ly, Z_EL))
    return Geometry(carbon=carbon, electrons=electrons, lx=lx, ly=ly, lz=LZ)


def atom_line(idx: int, atom_name: str, atom_type: str, x: float, y: float, z: float) -> str:
    return (
        f"HETATM {idx:5d} {atom_name:<5} "
        f"{'':<3} {'':<1} {atom_type:<5}"
        f"{x:10.5f}{y:10.5f}{z:10.5f}"
        f" {'':<5}{0:3d}{0:2d} {0.0:8.5f}\n"
    )


def write_geo(path: Path, geom: Geometry, label: str) -> None:
    nel = len(geom.electrons)
    ncarbon = len(geom.carbon)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("BIOGRF 200\n")
        handle.write(f"DESCRP {label}\n")
        handle.write("RUTYPE NORMAL RUN\n")
        handle.write("REMARK Jamil-style paired C+/e- graphene geometry\n")
        handle.write("REMARK fixed-charge MOLCHARGE branch (icharg=8)\n")
        handle.write(f"REMARK paired-electron radius={PAIR_RADIUS:.3f} A\n")
        handle.write(
            "FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)\n"
        )
        handle.write(
            f"CRYSTX {geom.lx:12.5f} {geom.ly:12.5f} {geom.lz:12.5f}"
            f" {90.0:12.5f} {90.0:12.5f} {90.0:12.5f}\n"
        )
        for idx, (x, y, z) in enumerate(geom.electrons, start=1):
            handle.write(atom_line(idx, "El", "El", x, y, z))
        for idx, (x, y, z) in enumerate(geom.carbon, start=nel + 1):
            handle.write(atom_line(idx, "C", "C", x, y, z))
        for i in range(1, nel + 1):
            handle.write(f"MOLCHARGE {i:4d} {i:4d}  -1.00\n")
        handle.write(f"MOLCHARGE {nel + 1:4d} {nel + ncarbon:4d}  {float(nel):6.2f}\n")
        handle.write("END\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--nx", type=int, default=7)
    parser.add_argument("--ny", type=int, default=4)
    parser.add_argument("--outdir", type=Path, default=Path("build/jamil_replication/geometries"))
    args = parser.parse_args()
    if args.nx <= 0 or args.ny <= 0:
        raise SystemExit("nx and ny must be positive")

    args.outdir.mkdir(parents=True, exist_ok=True)
    neutral = build_geometry(args.nx, args.ny, excess=False)
    excess = build_geometry(args.nx, args.ny, excess=True)
    ncarbon = len(neutral.carbon)
    write_geo(args.outdir / f"jamil_graphene_{ncarbon}c_{ncarbon}el.geo", neutral, f"jamil_graphene_{ncarbon}c_neutral")
    write_geo(args.outdir / f"jamil_graphene_{ncarbon}c_{ncarbon + 1}el.geo", excess, f"jamil_graphene_{ncarbon}c_excess")
    print(
        f"wrote ncarbon={ncarbon} neutral_el={len(neutral.electrons)} "
        f"excess_el={len(excess.electrons)} lx={neutral.lx:.5f} ly={neutral.ly:.5f}"
    )


if __name__ == "__main__":
    main()
