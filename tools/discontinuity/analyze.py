import sys, re, numpy as np
d = sys.argv[1]
rows = []
on = False
for l in open(d + "/stdout.txt"):
    if l.strip().startswith("Step"): on = True; continue
    if l.startswith("Loop time"): on = False
    if on:
        v = l.split()
        if len(v) == 19: rows.append([float(x) for x in v])
T = np.array(rows); n = len(T)
# dump
X, F = [], []
L = open(d + "/fd.dump").read().split("ITEM: TIMESTEP\n")[1:]
for b in L:
    ls = b.split("ITEM: ATOMS")[1].split("\n")[1:]
    a = np.array([[float(x) for x in s.split()] for s in ls if s.strip()])
    X.append(a[:, 2:5]); F.append(a[:, 5:8])
X = np.array(X); F = np.array(F)
pe = T[:, 1]; etot = T[:, 3]; terms = T[:, 4:]
print("steps", n, "etot span %.5f  start %.5f end %.5f" % (etot.max() - etot.min(), etot[0], etot[-1]))
# trapezoid work residual: dPE + 0.5 (F_t+F_t+1).dx   (F = -dPE/dx)
dx = X[1:] - X[:-1]
w = 0.5 * ((F[:-1] + F[1:]) * dx).sum(axis=(1, 2))
r = (pe[1:] - pe[:-1]) + w
print("work residual: median |r| %.2e  max |r| %.2e" % (np.median(abs(r)), abs(r).max()))
print("cum sum r end %.4f" % r.sum())
idx = np.argsort(-abs(r))[:15]
print("largest residual steps:", [(int(i), float("%.4g" % r[i])) for i in sorted(idx)])
np.save(d + "/r.npy", r)
# etot jumps
de = np.diff(etot)
print("etot step diffs: median |.| %.2e max %.2e; cumulative drift" % (np.median(abs(de)), abs(de).max()))
# term-wise: second difference of each term to see jumps
names = "eb ea elp emol ev epen ecoa ehb et eco ew ep efi eqeq eelect".split()
for k, nm in enumerate(names):
    t = terms[:, k]
    if np.all(t == 0): continue
    d2 = np.abs(t[2:] - 2 * t[1:-1] + t[:-2])
    print("%-6s max|d2| %.3e at step %d, median %.2e" % (nm, d2.max(), d2.argmax() + 1, np.median(d2)))
