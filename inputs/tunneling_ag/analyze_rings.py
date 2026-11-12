#!/usr/bin/env python3
"""A7: ring-polymer analysis for the tunneling demo.

Reads i-PI xyz bead trajectories <dir>/<prefix>.pos_<b>.xyz (b = 0..P-1; positions in Angstrom) for the TRANSFERRING electron
(an 'El' atom; default: the last El), projects onto the donor-acceptor axis and reports
  - histogram of bead positions along the axis (all beads, all frames)
  - histogram of centroid positions
  - fraction of bead density beyond the barrier midpoint (s > 0, the acceptor side when the electron starts on the donor side;
    sign is chosen so that the donor is at s < 0) and inside the forbidden region [-w, +w]
  - number of centroid crossings between the two sides (hysteresis +-h around s = 0)
  - radius of gyration per electron (mean, std over frames)
Axis: defined from two atoms of the FIRST frame (--axis-atoms i j, 0-based indices into the xyz atom list: donor then acceptor
reference atoms); r_mid = their midpoint, n_hat = (r_j - r_i)/|r_j - r_i|.

Usage: analyze_rings.py --axis-atoms 0 3 --out fig.png --title "m = 1 amu" classical=DIR1 pi=DIR2 [--discard 0.2] [--w 0.25] [--h 0.1]
Each label=DIR is overlaid in the figure; a JSON summary is written next to the figure.
"""
import argparse, json, pathlib, re, sys
import numpy as np


def read_xyz(path):
    """Return (symbols, frames[nframes, natoms, 3])."""
    syms, fr, cur = None, [], []
    lines = pathlib.Path(path).read_text().splitlines()
    i = 0
    while i < len(lines):
        n = int(lines[i].split()[0]); blk = lines[i + 2:i + 2 + n]
        if len(blk) < n:
            break  # truncated last frame
        w = [b.split() for b in blk]
        if syms is None: syms = [x[0] for x in w]
        fr.append([[float(v) for v in x[1:4]] for x in w]); i += n + 2
    return syms, np.array(fr)


def load_beads(d, prefix="sim"):
    d = pathlib.Path(d); files = sorted(d.glob(prefix + ".pos_*.xyz"), key=lambda p: int(re.findall(r"pos_(\d+)", p.name)[0]))
    if not files: raise FileNotFoundError("no %s.pos_*.xyz in %s" % (prefix, d))
    syms, first = read_xyz(files[0]); nf = len(first)
    beads = [first] + [read_xyz(f)[1][:nf] for f in files[1:]]
    nf = min(len(b) for b in beads)
    return syms, np.array([b[:nf] for b in beads])  # (P, nframes, natoms, 3)


def analyse(beads, syms, axis_atoms, el_index=None, discard=0.2, w=0.25, h=0.1, bins=80, srange=None, top=None):
    P, nf = beads.shape[:2]
    ref = beads[:, 0].mean(axis=0)                                  # first-frame centroids of every atom
    ri, rj = ref[axis_atoms[0]], ref[axis_atoms[1]]
    n = (rj - ri) / np.linalg.norm(rj - ri); rmid = 0.5 * (ri + rj)
    els = [k for k, s in enumerate(syms) if s == "El"]
    k = els[-1] if el_index is None else els[el_index]
    s_all = (beads[:, :, k, :] - rmid) @ n                          # (P, nf) bead coordinates along the axis
    c_all = beads[:, :, k, :].mean(axis=0)                           # centroid (nf, 3)
    s_c = (c_all - rmid) @ n
    rg = np.sqrt(((beads[:, :, k, :] - c_all[None]) ** 2).sum(-1).mean(axis=0))
    a = int(discard * nf); s_all, s_c, rg = s_all[:, a:], s_c[a:], rg[a:]
    half = 0.5 * np.linalg.norm(rj - ri)
    srange = srange or (-half - 1.0, half + 1.0)
    hb, edges = np.histogram(s_all.ravel(), bins=bins, range=srange, density=True)
    hc, _ = np.histogram(s_c, bins=bins, range=srange, density=True)
    # centroid crossings with hysteresis: state flips only when |s| > h on the other side
    state, cross = 0, 0
    for v in s_c:
        if state == 0:
            if abs(v) > h: state = 1 if v > 0 else -1
        elif state == -1 and v > h: state, cross = 1, cross + 1
        elif state == 1 and v < -h: state, cross = -1, cross + 1
    # statistical error: block averages (10 blocks) of the beyond-midpoint fraction
    def frac(mask_fn):
        b = np.array_split(np.arange(s_all.shape[1]), 10)
        v = np.array([mask_fn(s_all[:, i]).mean() for i in b if len(i)]); return float(v.mean()), float(v.std(ddof=1) / np.sqrt(len(v))) if len(v) > 1 else float("nan")
    fb, fbe = frac(lambda x: x > 0); ff, ffe = frac(lambda x: np.abs(x) < w)
    zones = {}
    if top is not None:   # zones split at the intrinsic barrier tops s = -top, +top (donor pit side s < -top, middle, far unit s > top)
        lab = lambda x: np.where(x < -top, 0, np.where(x > top, 2, 1))
        zb, zc = lab(s_all), lab(s_c)
        fz = lambda k: float((zb == k).mean())
        straddle = ((zb == 0).any(axis=0) & (zb > 0).any(axis=0)) | ((zb == 2).any(axis=0) & (zb < 2).any(axis=0))
        zones = dict(top=float(top), frac_beads_donor_side=fz(0), frac_beads_middle=fz(1), frac_beads_far_unit=fz(2),
                     frac_centroid_donor_side=float((zc == 0).mean()), frac_centroid_middle=float((zc == 1).mean()), frac_centroid_far_unit=float((zc == 2).mean()),
                     frac_frames_ring_straddling_a_top=float(straddle.mean()), centroid_zone_changes=int((np.diff(zc) != 0).sum()))
    return dict(zones=zones, P=int(P), nframes_used=int(s_all.shape[1]), axis=n.tolist(), r_mid=rmid.tolist(), half_gap=float(half),
                frac_beads_beyond_mid=fb, frac_beads_beyond_mid_err=fbe, frac_beads_forbidden=ff, frac_beads_forbidden_err=ffe,
                forbidden_halfwidth=w, centroid_frac_beyond_mid=float((s_c > 0).mean()), centroid_crossings=int(cross), hysteresis=h,
                rg_mean=float(rg.mean()), rg_std=float(rg.std()), bead_hist=hb.tolist(), centroid_hist=hc.tolist(), hist_edges=edges.tolist())


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("runs", nargs="+", help="label=DIR")
    ap.add_argument("--axis-atoms", nargs=2, type=int, required=True); ap.add_argument("--el-index", type=int, default=None)
    ap.add_argument("--out", required=True); ap.add_argument("--title", default=""); ap.add_argument("--prefix", default="sim")
    ap.add_argument("--top", type=float, default=None, help="barrier-top distance from the midpoint (A): adds zone fractions and straddling statistics")
    ap.add_argument("--discard", type=float, default=0.2); ap.add_argument("--w", type=float, default=0.25); ap.add_argument("--h", type=float, default=0.1)
    a = ap.parse_args()
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    res = {}
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6), sharey=False)
    for item in a.runs:
        lab, d = item.split("=", 1)
        syms, beads = load_beads(d, a.prefix)
        r = analyse(beads, syms, a.axis_atoms, a.el_index, a.discard, a.w, a.h, top=a.top); res[lab] = {k: v for k, v in r.items()}
        x = 0.5 * (np.array(r["hist_edges"][1:]) + np.array(r["hist_edges"][:-1]))
        if r["zones"]: print("             zones (donor|middle|far) beads %.3f|%.3f|%.3f  centroid %.3f|%.3f|%.3f  frames with a ring straddling a top %.3f  centroid zone changes %d" % (r["zones"]["frac_beads_donor_side"], r["zones"]["frac_beads_middle"], r["zones"]["frac_beads_far_unit"], r["zones"]["frac_centroid_donor_side"], r["zones"]["frac_centroid_middle"], r["zones"]["frac_centroid_far_unit"], r["zones"]["frac_frames_ring_straddling_a_top"], r["zones"]["centroid_zone_changes"]))
        ax[0].plot(x, r["bead_hist"], label="%s (P=%d)" % (lab, r["P"])); ax[1].plot(x, r["centroid_hist"], label=lab)
        print("%-12s P=%-3d beyond-mid %.4f+-%.4f  forbidden(|s|<%.2f) %.4f+-%.4f  centroid crossings %d  Rg %.3f+-%.3f A" % (
            lab, r["P"], r["frac_beads_beyond_mid"], r["frac_beads_beyond_mid_err"], a.w, r["frac_beads_forbidden"], r["frac_beads_forbidden_err"],
            r["centroid_crossings"], r["rg_mean"], r["rg_std"]))
    for k, t in ((0, "bead positions"), (1, "centroid positions")):
        ax[k].set_xlabel("s along donor-acceptor axis (A), barrier at 0"); ax[k].set_ylabel("density"); ax[k].set_title(t); ax[k].axvline(0, color="k", lw=0.5); ax[k].legend()
    fig.suptitle(a.title); fig.tight_layout(); fig.savefig(a.out, dpi=130)
    json.dump(res, open(pathlib.Path(a.out).with_suffix(".json"), "w"))


if __name__ == "__main__":
    main()
