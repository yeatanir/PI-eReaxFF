#!/usr/bin/env python3
"""Strict comparison of the exported (full-precision) ereac step-0 terms and gradients with the saved V2 oracle.

Usage: compare_phase0_full.py <ereac_export.out>
V2 reference: ereaxff_pimd_modular/build/jamil_replication/finite_difference_oracle_force/central/
Legacy 'force_work' of V2 relates to the gradient d (kcal/mol/A) by a constant; it is fitted and reported.
"""
import pathlib
import sys

import numpy as np

V2 = pathlib.Path(__file__).resolve().parents[3] / "ereaxff_pimd_modular/build/jamil_replication/finite_difference_oracle_force/central"
NAMES = ["bond", "overcoordination", "lone_pair", "molecule", "valence_penalty", "angle_conjugation", "hydrogen_bond",
         "torsion", "torsion_conjugation", "vdw", "coulomb", "charge", "external_field", "electron"]


def read_export(path):
    terms, grad, step = None, [], None
    for line in pathlib.Path(path).read_text().splitlines():
        w = line.split()
        if w[0] == "STEP":
            step = int(w[1])
            if step > 0:
                break
        elif step == 0 and w[0] == "TERMS":
            terms = np.array(w[1:], float)
        elif step == 0 and w[0] == "ESTRC":
            estrc = float(w[1])
        elif step == 0 and w[0] == "GRAD":
            grad.append([float(x) for x in w[2:5]])
    return terms, estrc, np.array(grad)


def main():
    terms, estrc, grad = read_export(sys.argv[1])
    v2 = {}
    for l in (V2 / "single_point_energy_terms.out").read_text().splitlines():
        if not l.startswith("#"):
            k, v = l.split()
            v2[k] = float(v)
    print("term                      ereac (kcal/mol)        V2 (kcal/mol)        abs diff")
    worst = 0.0
    for n, e in zip(NAMES, terms):
        d = abs(e - v2[n]); worst = max(worst, d)
        print(f"{n:24s} {e:20.12f} {v2[n]:20.12f} {d:12.3e}")
    print(f"{'total':24s} {estrc:20.12f} {v2['total']:20.12f} {abs(estrc - v2['total']):12.3e}")
    rows = [l.split() for l in (V2 / "single_point_forces.out").read_text().splitlines() if not l.startswith("#")]
    fw = np.array([[float(x) for x in r[1:4]] for r in rows])
    assert fw.shape == grad.shape, (fw.shape, grad.shape)
    c = -np.sum(fw * grad) / np.sum(grad * grad)   # force_work = -c * d
    resid = fw + c * grad
    print(f"\nforce_work = -c * gradient, fitted c = {c:.9e}  (4.184e-4 expected: kcal/mol/A -> amu A/fs^2)")
    print(f"max |force_work| = {np.abs(fw).max():.4e}; max residual = {np.abs(resid).max():.3e}; "
          f"relative max residual = {np.abs(resid).max() / np.abs(fw).max():.3e}")
    print(f"max |gradient| = {np.abs(grad).max():.6f} kcal/mol/A")


main()
