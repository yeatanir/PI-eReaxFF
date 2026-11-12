<p align="center">
  <img src="assets/logo.png" alt="PI-eReaxFF logo: a pi-shaped ring polymer of beads and springs inside an electron-density blob, next to the text eReaxFF" width="560">
</p>

# PI-eReaxFF

**What if an electron in a molecular-dynamics simulation could behave like a quantum particle, and still be cheap enough to simulate for nanoseconds?**

This project tries to make that real. It lets an electron tunnel through a barrier, spread out, and "leak" across a gap, inside a ReaxFF-style force field, using path-integral molecular dynamics.

> **Status in one line:** it works end to end for silver clusters with explicit electrons (Ag/El), and we have checked it against an independent quantum solver. It is not yet a general-purpose engine. [What works and what doesn't](#where-things-stand) is below. The repository is private and the license is undecided.

---

## The big picture

**ReaxFF** is a fast way to simulate chemistry: atoms, bonds that form and break, charges. It treats every atom classically.

**eReaxFF** adds one trick. Some electrons become particles of their own (type `El`) that move around between atoms. Now the model can describe things plain ReaxFF cannot: conduction electrons, electron transfer, an electron sitting in a gap.

But an electron is light, and a light particle is not a classical dot. It is smeared out, and it can slip through barriers it doesn't have the energy to climb. A normal simulation misses all of that.

**Path-integral molecular dynamics (PIMD)** fixes it with a neat trick. Replace each quantum particle with a **necklace of beads** joined by springs:

```
   classical particle            quantum particle (a ring polymer)

          •                            •───•───•
                                     •           •
     one point,                      •           •      the more spread out the necklace,
     one place                         •───•───•        the more "quantum" the particle

   a light particle = a big floppy necklace;  a heavy one = a tight little clump
```

Simulate the necklace with ordinary dynamics and you get the quantum statistics of the particle (zero-point motion, tunneling, delocalisation), for the cost of running many copies of the classical system.

**This repository connects the three pieces:**

```
        i-PI                     LAMMPS + Metal-Reax                eReaxFF force field
   (moves the beads,   <--->   (computes forces on every   <--->   (the physics: atoms, ions,
    keeps temperature)          bead, fast)                          explicit electrons)
```

- **i-PI** is the well-known PIMD driver. We added one feature: some atoms (the silver ions) can stay classical while the electron is a necklace.
- **Metal-Reax** is a fast ReaxFF engine written in C++. We added an eReaxFF mode, checked term by term against the original Fortran code.
- **The glue**: scripts, quantum solvers for cross-checking, and a pile of validation reports.

## Why bother?

The motivating problem is **electron leakage across thin barriers** (think silicon devices). Classical simulations say the electron stays put. Reality says it leaks. A force field with a PIMD electron should capture that. Before trusting it on silicon, we needed a clean test where the right answer is known. That is the demonstration below.

## What we have shown so far

**The test:** one electron between two silver clusters, 4 Å apart, with a narrow barrier in the middle. We computed its distribution two ways:

1. a **3-D Schrödinger solver** (the quantum answer, no sampling), and
2. **our PIMD**, running the real eReaxFF force field through i-PI and LAMMPS.

**The result:** they agree, for four different electron masses.

| Electron mass | How much more likely to be found *at the barrier* than a classical electron would be | Agreement with the solver |
|---|---|---|
| 4 amu (heavy) | about 3× (PIMD 2.8, solver 2.6) | within error |
| 1 amu | about 6× (PIMD 6.5, solver 5.9) | within error |
| 0.5 amu | about 9× (PIMD 8.9, solver 8.4) | needed 100 beads, 50 was not enough |
| 0.1 amu (very light) | about 11× (PIMD 11.0, solver 10.4) | within error; this is mostly zero-point smearing, not tunneling |

The lighter the electron, the more it spreads over the barrier, which is the expected quantum trend. A classical simulation sits at 1×.

Two honest notes from the experiment: the first barrier we tried (the force field's own) turned out to be so low that the electron simply isn't trapped, so we designed a better one with the solver first. And in these runs the other four electrons were allowed to move while the solver froze them, so a few percent of difference in the tails is expected.

## Along the way, things we found and fixed

- **A hidden jump in the original code**: a tiny energy jump whenever two silver atoms drift across a 5 Å cutoff, built into the original code. It made long runs with moving atoms drift. We found it, measured it, and added a "smooth" setting that removes it.
- **About a 2× speed-up** from removing a redundant step in our i-PI patch (P=32).
- **A 58 % speed-up** from one number (the socket polling delay), with bit-for-bit identical results.
- **A cost surprise**: one force evaluation takes about 10 microseconds. i-PI's communication overhead is 20× larger. Speed now lives in the plumbing, not the physics.

## Where things stand

| | |
|---|---|
| **Works and is validated** | Ag/El with fixed charges; eReaxFF terms match the original code to 1e-10 or better; i-PI classical-atoms patch passes its tests; PIMD agrees with the quantum solver (above) |
| **Works with caveats** | "Smooth" setting is for Ag/El only; periodic boxes need edges of at least 2 × bond cutoff + 1 Å |
| **Not done / known problems** | Fixed charges only (no charge equilibration with electrons); only C, El and Ag are supported; graphene runs collapse after about 110 fs (cause unknown); the i-PI patch is not yet submitted upstream; one test is at 63 of 64; scripts expect the original folder layout |
| **In progress** | A production-quality version of the tunneling study (more runs, error bars from independent runs) |

## Try it

Honest warning: the scripts were written for a particular folder layout and have not been re-run from a fresh clone, so expect to fix a few paths. The steps, in short:

1. Apply the i-PI patch (`ipi/patches/`) to i-PI at commit `93162d2`.
2. Build the pair style (`pair_style/`) as a LAMMPS plugin.
3. Run one electron between two silver clusters:
   `python3 inputs/tunneling_ag/run_barrier_lammps.py --help`
4. Run every check in one go: `python3 tools/run_all_checks.py`

Exact commands, flags and the full reproduction recipe are in **[docs/TECHNICAL.md](docs/TECHNICAL.md)**.

## Map of the repository

| Folder | What's in it |
|---|---|
| `pair_style/` | The Metal-Reax engine with the eReaxFF mode |
| `ipi/` | The i-PI patch (classical atoms next to quantum ones) |
| `drivers/`, `inputs/`, `tools/` | Glue code, run scripts, the Ag/El setup, solvers, one-command check runner |
| `validation/` | Every validation report, in order, with the evidence |
| `docs/` | Technical reference, the build plan, the audit of the earlier code, status log |
| `legacy_derived/` | Small patches against the original Fortran code. Restricted; see below |

## Good to know

- **License:** undecided. Some content derives from the original group code, so public release is not decided. Keep the repository private for now.
- **The force field** is the published Ag/El set from Evangelisti, Fichthorn and van Duin, *J. Chem. Phys.* **153**, 104106 (2020).
- Built on [i-PI](https://github.com/i-pi/i-pi), [LAMMPS](https://www.lammps.org) and [Metal-Reax](https://github.com/yeatanir/Metal-Reax).
