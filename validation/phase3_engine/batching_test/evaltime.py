import sys, time, numpy as np, os
HERE = os.getcwd()
sys.path.insert(0, "."); import lmpdrv
R = os.path.abspath("../.."); plug = R + "/forks/metal-reax-lammps/build-ns/plugin/reaxmetaladapterplugin.so"
pos0 = np.array([[float(v) for v in l.split()[1:4]] for l in open(HERE + "/gen/init.xyz").read().splitlines()[2:]])
rng = np.random.default_rng(1)
P = 72; beads = [pos0.copy() for _ in range(P)]
for b in beads: b[10] += rng.normal(0, 0.35, 3)
out = {}
for mode in ["run1", "run0"]:
    d = lmpdrv.LmpBatch(datadir=HERE + "/gen/client0", plugin=plug, v0=0.0, mode=mode)
    ps = [b / lmpdrv.BOHR for b in beads]
    for p in ps[:5]: d.one(np.eye(3), p)             # warm up
    d.t_lmp = 0; d.n = 0; t = time.perf_counter(); res = []
    for rep in range(10):
        res = [d.one(np.eye(3), p) for p in ps]
    wall = time.perf_counter() - t; out[mode] = res
    print(mode, "per eval: wall %.1f us, in-lammps(scatter+run+gather) %.1f us" % (wall / (10 * P) * 1e6, d.t_lmp / d.n * 1e6))
e1 = np.array([r[0] for r in out["run1"]]); e0 = np.array([r[0] for r in out["run0"]])
f1 = np.array([r[1] for r in out["run1"]]); f0 = np.array([r[1] for r in out["run0"]])
print("run1 vs run0: max |dE| %.3g Ha, max |dF| %.3g Ha/bohr; E range %.4g..%.4g" % (abs(e1 - e0).max(), abs(f1 - f0).max(), e0.min(), e0.max()))
