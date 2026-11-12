"""Baseline (socket) run of the reference setup under cProfile of the i-PI server. usage: prof_base.py P N [mode: base|var]"""
import sys, os, re, shutil, subprocess, time
sys.path.insert(0, os.getcwd()); import ab
P, N, kind = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
out = ab.HERE / ("prof_" + kind); out.mkdir(exist_ok=True)
if kind == "base":
    gen = ab.HERE / "gen"; shutil.copy(gen / "init.xyz", out / "init.xyz")
    if not (out / "client0").exists(): shutil.copytree(gen / "client0", out / "client0", symlinks=True)
    x = (gen / "input.xml").read_text(); x = re.sub(r"<checkpoint stride='\d+'", "<checkpoint stride='%d'" % N, x); x = re.sub(r"nbeads='\d+'", "nbeads='%d'" % P, x); x = x.replace("<batch_size> 72 </batch_size>", "<batch_size> %d </batch_size>" % P)
    x = re.sub(r"<total_steps>\d+</total_steps>", "<total_steps>%d</total_steps>" % N, x); (out / "input.xml").write_text(x)
    sock = re.search(r"<address>(bl\d+)</address>", x).group(1)
    srv = subprocess.Popen([sys.executable, str(ab.HERE / "prof_ipi.py")], cwd=out, env=ab.ENV, stdout=open(out / "ipi.log", "w"), stderr=subprocess.STDOUT)
    time.sleep(4)
    lm = subprocess.Popen([str(ab.ROOT / "external/install_lammps/bin/lmp"), "-in", "in.ipi", "-log", "none", "-nocite"], cwd=out / "client0", env=ab.ENV, stdout=open(out / "lmp.log", "w"), stderr=subprocess.STDOUT)
    bs = subprocess.Popen([sys.executable, str(ab.IPI / "bin/i-pi-py_driver"), "-m", "gaussian_barrier", "-u", "-a", sock + "b", "-o", "V0=%.15g,sigma=%.15g,rx=0,ry=0,rz=0,nx=0,ny=0,nz=1,first=10,last=10" % (ab.V0, ab.SIG)], cwd=out, env=ab.ENV, stdout=open(out / "drv.log", "w"), stderr=subprocess.STDOUT)
    t = time.time(); srv.wait(); wall = time.time() - t; lm.terminate(); bs.terminate()
else:
    ab.variant(P, 1, "prof_var"); sys.exit()
print("wall", wall)
