"""Donor-barrier-acceptor Ag/El systems for the tunneling demo (axis = z, Ag atoms on x=0 plane, non-periodic).

kind 'pair'  : Ag3 (L) + Ag3 (R) facing apex to apex, gap g between the two apex atoms, 5 electrons:
               L has 3 El (the apex El is the *transferring* electron, listed LAST so it is the last atom),
               R has 2 El on its two base atoms. Left end state (L neutral, R cation) and right end state
               (L cation, R neutral) are mirror images, so the two sites are degenerate by symmetry.
kind 'ion'   : Ag3 (L, 3 El) + one Ag+ (R, no El): asymmetric alternative.
Geo: MOLCHARGE lines BEFORE END (the legacy reader stops at END). El placed at the cation +0.001 A in x.
"""
import pathlib, shutil
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
PROJ = HERE.parents[1]
AG3 = np.array([[0.0, 0.0, 1.73645], [0.0, 1.41119, -0.86820], [0.0, -1.41119, -0.86820]])  # SI, first Ag3 listing
H_BASE = AG3[0, 2] - AG3[1, 2]  # apex-to-base height 2.60465 A


def build(kind, gap):
    """Return (symbols, positions[A], charges, ia_transfer (index of the transferring El), z_apexL, z_apexR)."""
    zl, zr = -gap / 2, +gap / 2
    left = np.array([[a[0], a[1], zl + (a[2] - AG3[0, 2])] for a in AG3])            # apex at zl, bases below
    ag = [left]
    if kind == "pair":
        right = np.array([[a[0], a[1], zr - (a[2] - AG3[0, 2])] for a in AG3])       # apex at zr, bases above
        ag.append(right)
    elif kind == "ion":
        ag.append(np.array([[0.0, 0.0, zr]]))
    else:
        raise ValueError(kind)
    ag = np.vstack(ag)
    n_ag = len(ag)
    el_idx = [1, 2]                                   # left base atoms
    if kind == "pair":
        el_idx += [4, 5]                              # right base atoms
    el_idx += [0]                                     # left apex El = transferring electron (last atom)
    el = np.array([ag[i] + np.array([0.001, 0.0, 0.0]) for i in el_idx])
    pos = np.vstack([ag, el])
    sym = ["Av"] * n_ag + ["El"] * len(el)
    q = np.array([1.0] * n_ag + [-1.0] * len(el))
    return sym, pos, q, len(sym) - 1, zl, zr


def write_case(path, kind, gap, shift=40.0):
    """Write a ereac case directory (control, ffield, fort.3/4/20/35, geo). Positions shifted by +shift in all axes."""
    path = pathlib.Path(path); path.mkdir(parents=True, exist_ok=True)
    sym, pos, q, ia, zl, zr = build(kind, gap)
    g = ["BIOGRF 200", "DESCRP tunneling_%s_gap%.2f" % (kind, gap), "RUTYPE NORMAL RUN", "REMARK donor-barrier-acceptor Ag/El",
         "FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)"]
    for i, (s, p) in enumerate(zip(sym, pos + shift), 1):
        g.append("HETATM %5d %-5s %-3s %1s %-5s%10.5f%10.5f%10.5f %-5s%3d%2d %8.5f" % (i, s, s, " ", s, *p, "", 0, 0, 0.0))
    for i, c in enumerate(q, 1):
        g.append("MOLCHARGE %4d %4d %6.2f" % (i, i, c))
    g.append("END")
    (path / "geo").write_text("\n".join(g) + "\n"); (path / "fort.3").write_text("\n".join(g) + "\n")
    ff = (PROJ / "inputs/ag_el/ffield_Ag_El").read_text()
    (path / "ffield").write_text(ff); (path / "fort.4").write_text(ff)
    ctl = (PROJ / "reference_data/phase0_jamil_225atom/np_export/control").read_text()
    (path / "control").write_text(ctl)
    (path / "fort.20").write_text("  0\n  0   0: Normal run\n"); (path / "fort.35").write_text("0\n")
    return sym, pos, q, ia, zl, zr
