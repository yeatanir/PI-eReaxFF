#!/usr/bin/env python3
"""Task 0.3: compare rebuilt-ereac energy terms (fort.73) with the saved V2 oracle terms.

Usage: python3 compare_phase0.py   (run from ereaxff_build/project)
"""
import pathlib

V2 = pathlib.Path(__file__).resolve().parents[3] / "ereaxff_pimd_modular/build/jamil_replication/finite_difference_oracle_force/central/single_point_energy_terms.out"
REF = pathlib.Path(__file__).resolve().parents[1] / "reference_data/phase0_jamil_225atom"
# fort.73 columns after the iteration number
COLS = ["Ebond", "Eatom", "Elp", "Emol", "Eval", "Ecoa", "Ehbo", "Etors", "Econj", "Evdw", "Ecoul", "Echarge", "Efield", "Eelect"]
# V2 term name -> fort.73 column
MAP = {"bond": "Ebond", "overcoordination": "Eatom", "lone_pair": "Elp", "molecule": "Emol", "valence_penalty": "Eval",
       "angle_conjugation": "Ecoa", "hydrogen_bond": "Ehbo", "torsion": "Etors", "torsion_conjugation": "Econj",
       "vdw": "Evdw", "coulomb": "Ecoul", "charge": "Echarge", "external_field": "Efield", "electron": "Eelect"}

def read_v2():
    d = {}
    for line in V2.read_text().splitlines():
        if line.startswith("#"):
            continue
        k, v = line.split()
        d[k] = float(v)
    return d

def read_fort73(path):
    last = [l for l in path.read_text().splitlines() if l.strip() and not l.strip().startswith("Iter")][-1].split()
    return dict(zip(COLS, map(float, last[1:15])))

v2 = read_v2()
cases = [("periodic, vpar38=0", "jamil_vpar38_0"), ("periodic, vpar38=10", "jamil_vpar38_10"),
         ("non-periodic, vpar38=0", "jamil_np_vpar38_0"), ("non-periodic, vpar38=10", "jamil_np_vpar38_10")]
res = {name: read_fort73(REF / d / "fort.73") for name, d in cases}
print("| term | V2 saved (kcal/mol) | " + " | ".join(n for n, _ in cases) + " |")
print("|---|---:|" + "---:|" * len(cases))
for term, col in MAP.items():
    row = [f"{res[n][col]:.2f}" for n, _ in cases]
    print(f"| {term} | {v2[term]:.2f} | " + " | ".join(row) + " |")
print()
for n, _ in cases:
    mx = max(abs(res[n][MAP[t]] - v2[t]) for t in MAP)
    print(f"max |ereac - V2| over the 14 terms, {n}: {mx:.3g} kcal/mol")
