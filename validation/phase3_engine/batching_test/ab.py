"""A/B timing: baseline (socket lmp client + batched python barrier driver, run_barrier_lammps.py) vs variant (i-PI ffdirect, one in-process LAMMPS for all beads + barrier).
usage: ab.py P N1 N2 REPS   -> wall(N2)-wall(N1) over (N2-N1) steps removes start-up. Interleaved baseline/variant."""
import os, re, subprocess, sys, time, statistics, shutil
from pathlib import Path
HERE = Path("/private/tmp/claude-501/-Users-yeatanir-PhD-eReaxFF-PIMD/db7a76bf-0f4f-46c7-b130-602710a4488b/scratchpad/bt")   # spaceless symlink to scratch/batch_test
ROOT = Path("/Users/yeatanir/PhD/eReaxFF PIMD/ereaxff_build"); IPI = ROOT / "forks/i-pi"; RUN = ROOT / "project/inputs/tunneling_ag/run_barrier_lammps.py"
ENV = dict(os.environ, OMP_NUM_THREADS="1", DYLD_LIBRARY_PATH=str(ROOT / "external/install_lammps/lib"), PYTHONPATH="%s:%s" % (ROOT / "external/lammps/python", IPI), PYTHONUNBUFFERED="1")
V0, SIG = 3 / 627.5094740631, 0.1 / 0.529177210903
def floor(P, N, tag): return variant(P, N, tag, dummy=True)
def uptime(): return os.getloadavg()[0]

def base(P, N, tag, extra=()):
    out = HERE / ("w_" + tag); t = time.time()
    r = subprocess.run([sys.executable, str(RUN), "--out", str(out), "--beads", str(P), "--steps", str(N), "--dt-fs", "0.125", "--el-mass", "1", "--start-z", "-0.3", "--seed", "2001",
        "--v0-kcal", "3", "--sigma-z", "0.1", "--pos-stride", "40", "--stride", "20", "--latency", "1e-5", "--timeout", "3000", *extra], env=ENV, capture_output=True, text=True)
    m = re.search(r"exit (\d+) after ([\d.]+) s", r.stdout); assert m and m.group(1) == "0", r.stdout[-500:] + r.stderr[-500:]
    return float(m.group(2)), out

def variant(P, N, tag, mode="run1", batch=None, extra_xml=None, dummy=False):
    out = HERE / ("w_" + tag); out.mkdir(exist_ok=True); gen = HERE / "gen"
    shutil.copy(gen / "init.xyz", out / "init.xyz")
    if not (out / "client0").exists(): shutil.copytree(gen / "client0", out / "client0", symlinks=True)
    x = (gen / "input.xml").read_text()
    x = re.sub(r"<ffsocket name='lmp'.*?</ffsocket>\s*", "", x, flags=re.S); x = re.sub(r"<ffsocket name='barrier'.*?</ffsocket>\s*", "", x, flags=re.S)
    ffd = ("<ffdirect name='lmp' threaded='False'><pes> custom </pes><pes_path>%s</pes_path><batch_size>%d</batch_size>"
           "<parameters>{datadir: %s, v0: %.16g, sigma: %.16g, atom: 10, mode: %s}</parameters></ffdirect>\n" % (HERE / "lmpdrv.py", batch or P, out / "client0", V0, SIG, mode))
    if dummy: ffd = "<ffdirect name='lmp' threaded='False'><pes> dummy </pes><batch_size>%d</batch_size></ffdirect>\n" % (batch or P)
    x = x.replace("<system>", ffd + "  <system>", 1).replace("<force forcefield='barrier'/>", "")
    x = re.sub(r"<checkpoint stride='\d+'", "<checkpoint stride='%d'" % N, x); x = re.sub(r"nbeads='\d+'", "nbeads='%d'" % P, x); x = re.sub(r"<total_steps>\d+</total_steps>", "<total_steps>%d</total_steps>" % N, x)
    if extra_xml: x = extra_xml(x)
    (out / "input.xml").write_text(x); t = time.time()
    r = subprocess.run([sys.executable, str(IPI / "bin/i-pi"), "input.xml"], cwd=out, env=ENV, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout[-800:] + r.stderr[-1500:]
    return time.time() - t, out

if __name__ == "__main__":
    P, N1, N2, REPS = map(int, sys.argv[1:5]); kinds = sys.argv[5].split(",") if len(sys.argv) > 5 else ["base", "var", "floor"]; res = {k: [] for k in kinds}
    for rep in range(REPS):
        for name, fn in [(k, {"base": base, "var": variant, "floor": floor}[k]) for k in kinds]:
            ta, _ = fn(P, N1, "%s_P%d_a" % (name, P)); tb, _ = fn(P, N2, "%s_P%d_b" % (name, P))
            sps = (N2 - N1) / (tb - ta); res[name].append(sps)
            print("rep %d %-4s P=%d: %.1f steps/s  (t_N1 %.1f s, t_N2 %.1f s, load1 %.1f)" % (rep, name, P, sps, ta, tb, uptime()), flush=True)
    print("MEDIAN P=%d: " % P + ", ".join("%s %.1f steps/s (%.2f ms/step)" % (k, statistics.median(v), 1e3 / statistics.median(v)) for k, v in res.items()))
