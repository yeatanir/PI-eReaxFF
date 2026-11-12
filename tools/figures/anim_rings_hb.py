"""Ring-polymer renderings for the HIGHER barrier (V0 = 8 kcal/mol, sigma = 0.10 A): same layout as anim_rings_ag.py (Ag atoms drawn), but
source runs are the Mac-Mini higher-barrier runs (scratch/hb_mini/runs, read-only, ring.npy = (beads, frames, xyz)). Imaginary-time necklaces,
NOT real-time motion. Writes anim_A_1amu_hb_ag.mp4, anim_A_0p5amu_hb_ag.mp4, anim_B_compare_hb_ag.mp4 (classical | 1 amu | 0.5 amu) and
anim_windows_hb.json. Frames are 25 / 18 / 50 fs apart (position stride 200 / 200 / 400), so windows are 4 ps. Run: python3 figures/anim_rings_hb.py"""
import sys, json, re
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import anim_rings as AR
import anim_rings_ag as AG
from anim_rings_ag import *

V0 = 8.0
AR.zV, AR.EV_bare, AR.EV = potential(v0=V0)       # Ring.update / setup read these module globals
HB = SCR / "hb_mini" / "runs"
RUNS_HB = {  # key: (run, label, window start search)  -- 1 amu seed 3001 (L start, 6 crossings), 0.5 amu seed 3011 (L start, 21 crossings)
    "cl": ("hb_cl_L_s3021", "classical, P=1"), 1.0: ("hb_m1_L_s3001", "1 amu, P=72"), 0.5: ("hb_m05_L_s3011", "0.5 amu, P=100")}
WIN = 4.0
INFO = {}
YMAX = 12
ZL = (-2.6, 2.6)


def load(key):
    d = HB / RUNS_HB[key][0]; xml = (d / "input.xml").read_text()
    dt = float(re.search(r"<timestep[^>]*>\s*([\d.eE+-]+)", xml).group(1)) * int(re.search(r"trajectory filename='pos' stride='(\d+)'", xml).group(1)) / 1000.0
    pos = np.load(d / "ring.npy").transpose(1, 0, 2).astype(float)      # (frames, beads, 3)
    return pos, dt


def window(key):
    pos, dt = load(key); n = int(round(WIN / dt)); zc = pos[:, :, 2].mean(1)
    if key == "cl":                                                      # no crossings: take the window with the closest approach to z = 0
        k = int(np.argmin(np.abs(zc[int(1 / dt):])) + int(1 / dt)); i0 = max(0, min(k - n // 2, len(zc) - n - 1)); w = (i0, i0 + n)
    else:
        w = pick_window(pos, dt, WIN, need_sides=0.12)
        if w is None: raise SystemExit("no crossing window for %s" % key)
    INFO[str(key)] = dict(run=RUNS_HB[key][0], P=int(pos.shape[1]), t0_ps=round(w[0] * dt, 3), t1_ps=round(w[1] * dt, 3), frame_dt_fs=round(dt * 1000, 1),
                          **{k: round(v, 3) for k, v in crossing_stats(pos, dt).items()})
    return swap(pos[w[0]:w[1] + 1]), dt


def colr(key):
    return CLASSICAL if key == "cl" else (MASS_COL[key] if key in MASS_COL else "#e69f00")


def _title(axt, key):
    axt.text(0.98, 0.9, "%s (imaginary-time necklace)" % RUNS_HB[key][1], transform=axt.transAxes, fontsize=6.5, va="top", ha="right", color="#555555")


def anim_A(key, path):
    fr, dt = window(key)
    fig = plt.figure(figsize=(7.2, 6.4)); gs = fig.add_gridspec(3, 1, height_ratios=[1, 1.25, 0.75], hspace=0.28, left=0.09, right=0.98, top=0.95, bottom=0.07)
    axt = fig.add_subplot(gs[0]); axb = fig.add_subplot(gs[1], sharex=axt); axo = fig.add_subplot(gs[2]); zoom_setup(axt, axb, ZL, (-1.0, 1.0), YMAX)
    r = Ring(axt, axb, colr(key), fr.shape[1]); r.beads.set_zorder(4); r.poly.set_zorder(3); pt = overview(axo, colr(key))
    tl = axt.text(0.02, 0.9, "", transform=axt.transAxes, fontsize=8, va="top"); _title(axt, key)
    axt.text(0.02, 0.72, "barrier 8 kcal/mol", transform=axt.transAxes, fontsize=6.5, color="#555555")
    def upd(i):
        r.update(fr, i); pt.set_data([fr[i][:, 2].mean()], [fr[i][:, 0].mean()]); tl.set_text("t = %.2f ps" % (i * dt)); return []
    FuncAnimation(fig, upd, frames=len(fr), interval=50).save(path, writer=FFMpegWriter(fps=20, bitrate=2000), dpi=140); plt.close(fig)


def anim_B(path):
    keys = ["cl", 1.0, 0.5]; data = {k: window(k) for k in keys}
    fig = plt.figure(figsize=(9.5, 4.8)); gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.5], hspace=0.08, wspace=0.28, left=0.06, right=0.99, top=0.9, bottom=0.12)
    rings = []; tls = []
    for c, k in enumerate(keys):
        axt = fig.add_subplot(gs[0, c]); axb = fig.add_subplot(gs[1, c], sharex=axt); zoom_setup(axt, axb, ZL, (-1.0, 1.0), YMAX, labels=(c == 0))
        rings.append(Ring(axt, axb, colr(k), data[k][0].shape[1])); rings[-1].beads.set_zorder(4); rings[-1].poly.set_zorder(3)
        axt.set_title(RUNS_HB[k][1], fontsize=8, color=colr(k) if k != "cl" else "k"); tls.append(axt.text(0.03, 0.9, "", transform=axt.transAxes, fontsize=7, va="top"))
        if c: axt.set_ylabel(""); axb.set_ylabel("")
    # frame spacing differs (25, 25, 50 fs): index each panel by physical time
    T = min(len(data[k][0]) * data[k][1] for k in keys); nfr = int(T / data["cl"][1]) if False else int(T / 0.025)
    def upd(i):
        t = i * 0.025
        for k, r, tl in zip(keys, rings, tls):
            fr, dt = data[k]; j = min(int(round(t / dt)), len(fr) - 1); r.update(fr, j); tl.set_text("t = %.2f ps" % (j * dt))
        return []
    FuncAnimation(fig, upd, frames=nfr, interval=50).save(path, writer=FFMpegWriter(fps=20, bitrate=2800), dpi=130); plt.close(fig)


if __name__ == "__main__":
    anim_A(1.0, HERE / "anim_A_1amu_hb_ag.mp4"); anim_A(0.5, HERE / "anim_A_0p5amu_hb_ag.mp4"); anim_B(HERE / "anim_B_compare_hb_ag.mp4")
    json.dump(INFO, open(HERE / "anim_windows_hb.json", "w"), indent=1); print(json.dumps(INFO, indent=1))
