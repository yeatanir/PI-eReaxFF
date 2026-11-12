# Does evaluating all beads of a step in one batch make the PIMD runs several times faster?

**Verdict (measured).** Yes, about 5-6x at the reference setup, and the prototype is already close to its ceiling. Removing the per-bead socket round trip (one in-process LAMMPS inside i-PI, `ffdirect`, all P beads in one call, Gaussian barrier included) took the 1 amu setup (P=72, dt 0.125 fs, `s3_m1_d1_mob_*` settings) from **50.4 to 256 steps/s (5.1x)**, P=100 from 35.7 to 189 (5.3x), P=220 from 14.8 to 91 (6.2x). The forces are identical to 2.5e-12 kcal/mol/A on 793 bead configurations, and a same-seed 400-step trajectory agrees to 3e-6 eV. The README statement "several times more speed" holds. The README statement that the force field is only ~4 % of the cost needs a small correction (below).

Labels: **[measured]** this test; **[derived]** from measurements; **[estimate]** not measured.

## Setup [measured]
- Machine: MacBook, load average 11-17 the whole time (5 Stage 3 runs going), so absolute steps/s are lower than on an idle machine; baseline and variants were run **interleaved, 3 repeats**, medians below, spread 5.0-5.3x (P=72). My test used 1 core for the variant and the normal 3 processes for the baseline.
- Baseline = the reference run recipe exactly (`inputs/tunneling_ag/run_barrier_lammps.py`, `--beads P --dt-fs 0.125 --el-mass 1 --start-z -0.3 --seed 2001 --v0-kcal 3 --sigma-z 0.1 --latency 1e-5`, Ag fixed, 4 classical El): i-PI server + one LAMMPS client over a unix socket (one request per bead) + a Python Gaussian-barrier driver on a second batched socket.
- Variant = same input.xml with both `ffsocket`s replaced by one `ffdirect` (`pes=custom`, `lmpdrv.py`): the python `lammps` module holding one LAMMPS instance (same `lmp` build, same Metal-Reax plugin), for each bead scatter positions, `run 1 pre no post no`, gather forces and energy, add the barrier (same formula and units as `gaussian_barrier.py`). No i-PI source was touched.
- Steps/s = (N2-N1)/(t(N2)-t(N1)) with N1,N2 = 1000,3000 (P=72), 600,2000 (P=100), 300,900 (P=220), which removes start-up. "Floor" = same i-PI dynamics with a zero-force `dummy` ffdirect, no sockets and no LAMMPS: the cost of i-PI itself.

## Timings [measured, medians of 3 interleaved repeats]
| P | baseline (socket) | in-process batch | i-PI floor (zero force) | speedup batch/baseline | bound floor/baseline |
|---|---|---|---|---|---|
| 72 | 50.4 steps/s (19.8 ms/step) | 256.2 (3.90 ms) | 464.3 (2.15 ms) | **5.1x** | 9.2x |
| 100 | 35.7 (28.0 ms) | 188.6 (5.30 ms) | 346.2 (2.89 ms) | **5.3x** | 9.7x |
| 220 | 14.8 (67.5 ms) | 91.2 (10.96 ms) | 144.0 (6.94 ms) | **6.2x** | 9.7x |

## Where the 19.8 ms/step at P=72 goes [derived from the three runs]
| Part | ms/step | share |
|---|---|---|
| i-PI Python itself (propagation, PILE-L thermostat, bookkeeping), floor run | 2.15 | 11 % |
| LAMMPS force + barrier, in process (3.90 - 2.15) = 24 us per bead incl. ctypes scatter/gather | 1.75 | 9 % |
| Everything the sockets add (two socket round trips per bead, polling, serialising, the extra Python barrier-driver process, waits) | 15.9 | **80 %** |
Direct in-process timing of the force call alone (`evaltime.py`, P=72 random ring configurations): 22 us per bead with scatter/gather, 33 us with a full `run 0`. The "about 10 us" in the README is the bare force and is consistent: the extra ~12 us is the Python module's data copy. So the physics is **about 4 % (bare force) to 9 % (as called from Python)** of the baseline cost; the README's "about 4 %" is the low end, state it as "4 to 9 %".
A cProfile of the baseline server (1500 steps) agrees qualitatively: the largest entries are `time.sleep`, `select`, `socket send/recv` and the thread pool that distributes requests, with the PILE thermostat step (about 13 % of profile time) the biggest pure-numerics item. (A first profile was polluted by a per-step checkpoint in my generated input and was discarded; every number in this note uses a checkpoint only at the end, like the production runs.)

## Upper bound [derived]
T_now / (T_floor + P x T_force): with T_force = 10.4 us (bare) the bound at P=72 is 19.8 / (2.15 + 0.75) = **6.8x**; with the 24 us per bead the prototype actually pays it is 5.1x, which is what was measured. So a faster in-process call (C-level, no ctypes copies) could add at most ~30 % more; beyond that only i-PI's own Python cost (11 % of baseline) limits things, up to ~9x.
Cores: the baseline uses about 1.4 cores at 50 steps/s (a `ps` snapshot of Stage 3 runs: i-PI ~120 %, lmp ~16 %, barrier driver ~4 %); the variant is one process. At equal core count the gain is therefore larger than the per-run ratio, roughly 7x per core **[estimate]**.

## Correctness [measured]
1. Forces: 793 bead configurations taken from a baseline run (`--dump-forces`, P=72, 10 steps): max |F_in-process - F_socket client| = **2.5e-12 kcal/mol/A** (force scale up to 51), for both `run 1` and `run 0` modes (`corrcheck.py`). `run 1 pre no post no` and `run 0` give identical energies and forces on 72 random rings (max diff 0).
2. Same-seed trajectory, 400 steps, P=72, i-PI outputs every 20 steps: |diff| conserved 2.9e-6 eV (of 16.14), potential 3.1e-6 eV (of 34.9), kinetic_md 2.5e-7, kinetic_cv 4e-9, temperature 2e-6 K. Not bitwise: i-PI's unit constants for the socket client differ from mine in the 7th digit. Chaotic divergence over long runs is expected, so equivalence is statistical, as between any two seeds.

## What production use would still need [not done]
- A `--inproc` option in `run_barrier_lammps.py` that writes the `ffdirect` input (currently hand-generated in `ab.py`) and keeps the ring.npy reduction, run_info and `CONTROL_USED` logging unchanged.
- `plugin load` cannot take a path with spaces: the test symlinks the plugin as `plug.so` inside the data dir. A production wrapper should do the same.
- Restart: the driver is stateless (positions come from i-PI each call), so i-PI `RESTART` should work, but I did not test a restart. Thermostat (PILE-L), `fixatoms`, `classical_atoms` and the observables were exercised in the 400-step comparison.
- A long statistical comparison (density at the barrier, 50 ps) of in-process against socket runs before using it for results; I only compared 400 steps.
- Multi-run: one process per run, no socket files. i-PI-side overheads (11 % floor, thermostat) become the limit; moving them to NumPy-vectorised or compiled code is the next step after batching, **[estimate]** up to ~9x total.

Files: `batching_test/lmpdrv.py` (the PES), `ab.py` (A/B harness), `evaltime.py`, `corrcheck.py`, `prof_base.py`, `prof_ipi.py`, `ab2.log` (P=72 floor, P=100, P=220 timings). Scratch work: `scratch/batch_test/` (paths in the scripts point there).
