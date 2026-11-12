#!/usr/bin/env python3
"""Package selected frames of an ereac export (STEP/TERMS/ESTRC/GRAD/POS blocks) into one JSON reference file.
Usage: make_frames.py <ereac_export.out> <out.json> <step0:step1:stride> [<step0:step1:stride> ...]"""
import json, sys
exp, out = sys.argv[1:3]
sel = set()
for spec in sys.argv[3:]:
    a, b, s = map(int, spec.split(":")); sel.update(range(a, b + 1, s))
frames, cur = [], None
for l in open(exp):
    w = l.split()
    if w[0] == "STEP":
        cur = {"step": int(w[1]), "grad": [], "pos": []}; frames.append(cur)
    elif cur is None: continue
    elif w[0] == "TERMS": cur["terms"] = [float(x) for x in w[1:15]]
    elif w[0] == "ESTRC": cur["estrc"] = float(w[1])
    elif w[0] == "GRAD": cur["grad"].append([float(x) for x in w[2:5]])
    elif w[0] == "POS": cur["pos"].append([float(x) for x in w[2:5]])
frames = [f for f in frames if f["step"] in sel and f["pos"]]
json.dump({"terms_order": ["bond","overcoordination+undercoordination","lone_pair","molecule","valence+penalty","angle_conjugation","hydrogen_bond","torsion","torsion_conjugation","vdw","coulomb","charge","external_field","electron"],
           "units": {"energy": "kcal/mol", "gradient": "kcal/mol/A (force = -gradient)", "positions": "A"},
           "frames": frames}, open(out, "w"))
print(len(frames), "frames written")
