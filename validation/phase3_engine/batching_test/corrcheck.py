import sys, os, numpy as np
HERE = os.getcwd(); sys.path.insert(0, HERE); import lmpdrv
L = open(HERE + "/corr_base/client0/fd.dump").read().split("\n"); fr = []; i = 0
while i < len(L) and L[i].startswith("ITEM: TIMESTEP"):
    n = int(L[i + 3]); rows = [list(map(float, l.split())) for l in L[i + 9:i + 9 + n]]; fr.append(np.array(rows)); i += 9 + n
print("frames", len(fr))
for mode in ["run1", "run0"]:
    d = lmpdrv.LmpBatch(datadir=HERE + "/gen/client0", mode=mode); dE = dF = 0
    for a in fr:
        e, f = d.lammps_eval(a[:, 1:4]); dF = max(dF, abs(f - a[:, 4:7]).max())
    print(mode, "max |F_inproc - F_socket-client dump| = %.3g kcal/mol/A (force scale up to %.1f)" % (dF, max(abs(a[:, 4:7]).max() for a in fr)))
