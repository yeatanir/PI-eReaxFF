#!/usr/bin/env python3
"""Conserved-quantity bookkeeping of the Stage 2 runs (i-PI sim.out column 'conserved{electronvolt}', every 20 steps, PILE-L + fixed Ag + 4 classical El + one ring).
usage: stage2_conserved.py <dir> [<dir> ...]  -> per run: length, std and span (kT) after 1 ps, linear drift (kT per ns, +- from 2 halves), temperature mean."""
import sys, numpy as np
KT = 8.617333e-5 * 300.0
for d in sys.argv[1:]:
    a = np.loadtxt(d + "/sim.out"); t = a[:, 1] / 1000.0; c = a[:, 2]; T = a[:, 6]; m = t >= 1.0; t, c, T = t[m], c[m], T[m]
    sl = lambda x, y: np.polyfit(x, y, 1)[0]; h = len(t) // 2
    print("%-34s %6.1f ps  std %.3g kT  span %.3g kT  drift %+.3g kT/ns (halves %+.3g / %+.3g)  <T>_md %.1f K" % (d.split('/')[-1], t[-1], c.std() / KT, (c.max() - c.min()) / KT,
          sl(t, c) / KT, sl(t[:h], c[:h]) / KT, sl(t[h:], c[h:]) / KT, T.mean()))
