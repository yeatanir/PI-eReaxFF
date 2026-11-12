"""i-PI Python-driver PES whose forces come from the exact rebuilt legacy `ereac` (single point per call, via files).

Env: EREAC_BIN (export-build binary), EREAC_CASE (dir with control, ffield, geo, fort.20, fort.35 of the system; positions in geo are replaced),
     EREAC_WORKERS (parallel calls for batched requests, default 3), EREAC_TMP (scratch root).
Slow (process + files per call) but it IS the legacy code: use for validating the PIMD chain on small systems.
Positions are written with f10.5 (legacy geo format): 1e-5 A quantisation limits energy conservation to ~1e-3 kcal/mol.
Units: i-PI bohr/Hartree <-> A, kcal/mol.
"""
import os, re, shutil, subprocess, tempfile, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
from ipi.pes.dummy import Dummy_driver

__DRIVER_NAME__ = "ereac"
__DRIVER_CLASS__ = "EReac_driver"

BOHR_A = 0.529177210903
KCAL_PER_HARTREE = 627.5094740631  # exact-ish; ereac energies are kcal/mol


def _setkey(t, key, val):
    p = re.compile(r"^\s*\S+(\s+%s\b.*)$" % key, re.M)
    if p.search(t):
        return p.sub(lambda m: "%7s%s" % (val, m.group(1)), t, count=1)
    return t.rstrip("\n") + "\n%7s %s\n" % (val, key)


class EReac_driver(Dummy_driver):
    def __init__(self, *args, **kwargs):
        self.bin = Path(os.environ["EREAC_BIN"]).resolve()
        self.case = Path(os.environ["EREAC_CASE"]).resolve()
        self.workers = int(os.environ.get("EREAC_WORKERS", "3"))
        self.tmp = Path(os.environ.get("EREAC_TMP", tempfile.gettempdir())) / ("ereac_drv_%d" % os.getpid())
        self.tmp.mkdir(parents=True, exist_ok=True)
        self.geo = (self.case / "geo").read_text().split("\n")
        self.nat = sum(1 for l in self.geo if l.startswith("HETATM"))
        ctl = (self.case / "control").read_text()
        for k, v in (("imdmet", "3"), ("iravel", "0"), ("nmdit", "1"), ("iout1", "1"), ("imetho", "0")):
            ctl = _setkey(ctl, k, v)
        self.control = ctl
        self.pool = ThreadPoolExecutor(self.workers)
        self.ncall = 0
        super().__init__(*args, **kwargs)

    def _one(self, idx, pos_bohr):
        d = self.tmp / ("w%d" % idx)  # one directory per task: tasks run concurrently
        shutil.rmtree(d, ignore_errors=True)  # ereac leaves restart/state files; every call must start clean
        d.mkdir(exist_ok=True)
        x = np.asarray(pos_bohr).reshape(-1, 3) * BOHR_A
        out, i = [], 0
        for l in self.geo:
            if l.startswith("HETATM"):
                l = l[:30] + "".join("%10.5f" % v for v in x[i]) + l[60:]
                i += 1
            out.append(l)
        geo = "\n".join(out)
        (d / "geo").write_text(geo); (d / "fort.3").write_text(geo)
        ff = (self.case / "ffield").read_text()
        (d / "ffield").write_text(ff); (d / "fort.4").write_text(ff)
        shutil.copy(self.case / "fort.20", d / "fort.20"); shutil.copy(self.case / "fort.35", d / "fort.35")
        (d / "control").write_text(self.control)
        for attempt in range(80):
            (d / "ereac_export.out").unlink(missing_ok=True)
            r = subprocess.run([str(self.bin)], cwd=d, capture_output=True, text=True, timeout=300)
            if "dyld" in r.stderr:
                time.sleep(0.2 + 0.05 * attempt); continue
            break
        else:
            raise RuntimeError("ereac could not be started (dyld)")
        e, g, step = None, [], None
        for l in (d / "ereac_export.out").read_text().splitlines():
            w = l.split()
            if w[0] == "STEP":
                step = int(w[1])
            elif step == 0 and w[0] == "ESTRC":
                e = float(w[1])
            elif step == 0 and w[0] == "GRAD":
                g.append([float(v) for v in w[2:5]])
        g = np.array(g)
        if e is None or g.shape != (self.nat, 3) or not np.all(np.isfinite(g)):
            raise RuntimeError("bad ereac result: " + r.stdout[-300:])
        # kcal/mol, kcal/mol/A -> Hartree, Hartree/bohr ; force = -gradient
        return e / KCAL_PER_HARTREE, (-g * BOHR_A / KCAL_PER_HARTREE)

    def compute_structure(self, cell, pos):
        e, f = self._one(0, pos)
        return e, f.reshape(np.asarray(pos).shape), cell * 0.0, "nada"

    def compute(self, cell, pos):
        if not isinstance(cell, list):
            return self.compute_structure(cell, pos)
        res = list(self.pool.map(lambda a: self._one(*a), list(enumerate(pos))))
        return [(e, f.reshape(np.asarray(p).shape), np.asarray(c) * 0.0, "nada") for (e, f), p, c in zip(res, pos, cell)]
