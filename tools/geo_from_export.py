#!/usr/bin/env python3
"""Write a new geo (same HETATM format) with positions taken from the LAST 'POS' block of an ereac export file.
Usage: geo_from_export.py <ereac_export.out> <template_geo> <out_geo> [frame_index]   (frame_index counts POS blocks, default last)"""
import sys
exp, tmpl, out = sys.argv[1:4]
frame = int(sys.argv[4]) if len(sys.argv) > 4 else -1
blocks, cur = [], None
prev = None
for l in open(exp):
    w = l.split()
    if w[0] == "POS":
        if prev != "POS":
            cur = []; blocks.append(cur)
        cur.append([float(x) for x in w[2:5]])
    prev = w[0]
pos = blocks[frame]
lines, i = [], 0
for l in open(tmpl).read().split("\n"):
    if l.startswith("HETATM"):
        l = l[:30] + "".join("%10.5f" % x for x in pos[i]) + l[60:]
        i += 1
    lines.append(l)
assert i == len(pos), (i, len(pos))
open(out, "w").write("\n".join(lines))
