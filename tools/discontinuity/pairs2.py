import sys, numpy as np
d = sys.argv[1]; thr = float(sys.argv[2]) if len(sys.argv) > 2 else 0.008
r = np.load(d + "/r.npy")
L = open(d + "/fd.dump").read().split("ITEM: TIMESTEP\n")[1:]
X = []; typ = None
for b in L:
    ls = b.split("ITEM: ATOMS")[1].split("\n")[1:]
    a = np.array([[float(x) for x in s.split()] for s in ls if s.strip()])
    X.append(a[:, 2:5]); typ = a[:, 1].astype(int)
X = np.array(X); n = X.shape[1]
nm = {1: "Ag", 2: "El"}
D = np.linalg.norm(X[:, :, None] - X[:, None], axis=-1)
iu = np.triu_indices(n, 1)
ev = np.where(abs(r) > thr)[0]
print("events", len(ev), "of", len(r), "thr", thr, " sign + :", (r[ev] > 0).sum(), " - :", (r[ev] < 0).sum())
print("event |r| values: min %.4f max %.4f mean %.4f std %.4f" % (abs(r[ev]).min(), abs(r[ev]).max(), abs(r[ev]).mean(), abs(r[ev]).std()))
for i in ev[:40]:
    # pairs whose distance changed most across the step (candidates) -- report ones in [3,7] straddling a round value
    out = []
    for a, b in zip(*iu):
        d0, d1 = D[i, a, b], D[i + 1, a, b]
        out.append((abs(d1 - d0), a, b, d0, d1))
    print("step %d r=%+.4f" % (i, r[i]), end=" | ")
    # list pairs that straddle 5.7, 5.0, 4.0, 3.0 etc.? print which of the cutoff candidates are straddled
    st = []
    for a, b in zip(*iu):
        d0, d1 = D[i, a, b], D[i + 1, a, b]
        for c in np.arange(1.5, 8.0, 0.05).round(2):
            if (d0 - c) * (d1 - c) < 0: st.append("%s%d-%s%d@%.1f(%.3f->%.3f)" % (nm[typ[a]], a, nm[typ[b]], b, c, d0, d1))
    print(" ".join(st[:6]))
