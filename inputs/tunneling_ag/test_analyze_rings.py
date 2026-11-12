"""Synthetic checks of analyze_rings.py with known answers (run: python3 -m pytest -q test_analyze_rings.py)."""
import pathlib, sys
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import analyze_rings as ar
from math import erf, sqrt


def write(tmp, beads, syms):
    P, nf = beads.shape[:2]
    for b in range(P):
        out = []
        for f in range(nf):
            out += [str(len(syms)), "# synthetic"] + ["%s %.6f %.6f %.6f" % (s, *beads[b, f, i]) for i, s in enumerate(syms)]
        (tmp / ("sim.pos_%d.xyz" % b)).write_text("\n".join(out) + "\n")


def synth(P=8, nf=6000, sig=0.5, c0=-1.0, seed=0, hops=None):
    rng = np.random.default_rng(seed)
    syms = ["Av", "Av", "El"]                       # donor Ag at z=-2, acceptor Ag at z=+2, electron
    beads = np.zeros((P, nf, 3, 3)); beads[:, :, 0, 2] = -2.0; beads[:, :, 1, 2] = 2.0
    cen = np.full(nf, c0, float)
    if hops:
        cen = np.where(((np.arange(nf) // hops) % 2) == 0, c0, -c0)
    off = rng.normal(0, sig, (P, nf, 3)); off -= off.mean(axis=0)   # ring offsets, zero centroid
    beads[:, :, 2, :] = off; beads[:, :, 2, 2] += cen[None]
    return syms, beads


def test_known_fractions_rg_and_zero_crossings(tmp_path):
    P, sig = 8, 0.5
    syms, beads = synth(P=P, sig=sig, c0=-1.0)
    write(tmp_path, beads, syms)
    s, b = ar.load_beads(tmp_path)
    r = ar.analyse(b, s, (0, 1), discard=0.0, w=0.25, h=0.1)
    sd = sig * sqrt((P - 1) / P)                          # std of a recentred bead coordinate
    exp_beyond = 0.5 * (1 + erf((0 - (-1.0)) / (sd * sqrt(2)))) ; exp_beyond = 1 - exp_beyond
    exp_forb = 0.5 * (erf((0.25 + 1) / (sd * sqrt(2))) - erf((-0.25 + 1) / (sd * sqrt(2))))
    assert abs(r["frac_beads_beyond_mid"] - exp_beyond) < 4 * max(r["frac_beads_beyond_mid_err"], 2e-3), (r["frac_beads_beyond_mid"], exp_beyond)
    assert abs(r["frac_beads_forbidden"] - exp_forb) < 4 * max(r["frac_beads_forbidden_err"], 2e-3)
    assert r["centroid_crossings"] == 0
    assert abs(r["rg_mean"] - sig * sqrt(3 * (P - 1) / P)) < 0.01                       # <Rg^2>^(1/2) (mean slightly lower: <sqrt> < sqrt<>)
    assert abs(sum(np.array(r["bead_hist"]) * np.diff(r["hist_edges"])) - 1.0) < 1e-6


def test_crossings_with_hysteresis(tmp_path):
    syms, beads = synth(P=4, nf=2000, sig=0.05, c0=-1.0, hops=250)           # 8 blocks -> 7 sign changes of the centroid
    write(tmp_path, beads, syms)
    s, b = ar.load_beads(tmp_path)
    r = ar.analyse(b, s, (0, 1), discard=0.0, h=0.1)
    assert r["centroid_crossings"] == 7, r["centroid_crossings"]


def test_cli_and_figure(tmp_path):
    import subprocess
    syms, beads = synth(P=4, nf=500)
    d = tmp_path / "run"; d.mkdir(); write(d, beads, syms)
    out = tmp_path / "f.png"
    p = subprocess.run([sys.executable, str(pathlib.Path(ar.__file__)), "--axis-atoms", "0", "1", "--out", str(out), "pi=" + str(d)], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    assert out.exists() and out.with_suffix(".json").exists()


def test_zone_fractions_and_straddling(tmp_path):
    P, sig = 8, 0.5
    syms, beads = synth(P=P, sig=sig, c0=-1.0)                       # donor at z=-2, acceptor z=+2, ring centred at s=-1
    write(tmp_path, beads, syms)
    s, b = ar.load_beads(tmp_path)
    r = ar.analyse(b, s, (0, 1), discard=0.0, top=1.5)               # tops at s = -1.5 and +1.5
    sd = sig * sqrt((P - 1) / P)
    Phi = lambda x: 0.5 * (1 + erf(x / sqrt(2)))
    exp = (Phi((-1.5 + 1.0) / sd), Phi((1.5 + 1.0) / sd) - Phi((-1.5 + 1.0) / sd), 1 - Phi((1.5 + 1.0) / sd))
    z = r["zones"]
    got = (z["frac_beads_donor_side"], z["frac_beads_middle"], z["frac_beads_far_unit"])
    assert all(abs(a - e) < 0.01 for a, e in zip(got, exp)), (got, exp)
    assert abs(sum(got) - 1.0) < 1e-12 and z["frac_frames_ring_straddling_a_top"] > 0.5
