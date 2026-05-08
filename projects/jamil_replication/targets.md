# Jamil Paper Replication Targets

Paper: Md Jamil Hossain, Gorakh Pawar, Adri C. T. van Duin,
"Development and Applications of an eReaxFF Force Field for Graphitic Anodes of
Lithium-Ion Batteries", Journal of The Electrochemical Society 169, 110540
(2022).

Local PDF: `/Users/yeatanir/PhD/PIMD/Jamil_paper.pdf`

## Method Baseline To Reproduce First

The paper's simulations are classical eReaxFF, not PIMD.

| Quantity | Paper Value |
|---|---:|
| electron representation | C+/e- quasi-Drude-like pairs |
| pseudo-electron mass | 1.0 amu |
| timestep | 0.1 fs |
| thermostat | Berendsen |
| Berendsen damping | 100 fs |
| ensemble | NVT |
| initial relaxation | 1 K NVT relaxation |
| production temperatures | 300 K, 400 K, 500 K |
| charge model | fixed C+/e- pairs can avoid ACKS2 in their ADF setup |

Important divergence in our current driver:

- Current eReaxFF adapter still goes through the legacy oracle interface and
  force export path.
- Our current paired graphene baseline is 60 C + 60 El, not the exact 112 C
  paper cell.
- The mixed BAOAB driver currently treats electron groups as PIMD-capable; a
  `P=1` run is the near-classical stand-in until a true classical electron group
  path is added.

## Target Matrix

| Priority | Paper Result | System | Conditions | Observable | Local First Step |
|---:|---|---|---|---|---|
| 1 | Fig. 3, in-plane electron mobility | single-layer graphene, one e- per C, plus one excess e- | 300 K, 50 ps | excess electron host/path visits many C atoms | reproduce qualitatively on 60 C + 61 El first |
| 2 | Fig. 4, out-of-plane graphite restriction | four graphene layers, 32 C/layer, one excess e- in layer 2 | 300 K, 205 ps | excess e- stays in source layer | generate 4-layer 128 C + 129 El geometry |
| 3 | Fig. 4, temperature effect | same graphite stack | 400 K and 500 K | layer-transfer events increase | short high-T stress tests |
| 4 | Table II, voltage/excess loading | graphite, layer 2 loaded with 2-5 excess e- | 300 K, 200 ps | count out-of-plane electron events | add multiple excess electrons after target 2 |
| 5 | Fig. 5-7, sp3 linker | bilayer graphene connected by sp3 C-C linker | 300 K, 2/3/5 excess e- | interlayer transfer through linker vs vacuum | later geometry generator |
| 6 | Fig. 8-9, vacancy localization | graphene with single C vacancy, one excess e- | 300 K, 50 ps | electron localizes near vacancy/neighbors; paper reports 28% near defect zone | second single-layer target after mobility |
| 7 | Fig. 10, Li+ electron capture | graphene with one excess e- and Li+ | 300 K, 800 ps | Li+ captures excess e- near surface at about 3 ps | later, requires Li setup validation |
| 8 | Fig. 11-12, Li plating | graphene plus 5 Li+ and 5 excess e- | 300 K | Li neutralization and Li-Li clustering/plating | later, after Li target 7 |

## First Replication Target

Target: Fig. 3-style single-layer graphene in-plane excess-electron mobility.

Start smaller than the paper:

- use existing Track B2 paired graphene geometry as baseline;
- make a 60 C + 61 El charged variant;
- use `P = 1`, electron mass `1.0 amu`, `dt = 0.1 fs`, `T = 300 K`;
- use Berendsen `tau = 100 fs` for carbon and electrons if the driver supports
  true classical electrons, otherwise use a documented `P=1` PIMD stand-in;
- keep guard rails enabled for short tests.

Required analysis outputs:

- nearest-host carbon ID for the extra electron per frame;
- unique host count;
- hop count;
- nearest-host residence histogram;
- path-density XYZ/CSV for visualization.

Initial acceptance:

- 0-step and 1-step runs initialize with finite forces;
- generated oracle geos include `CRYSTX`;
- no near-zero El-El or El-C distances;
- timing summary is recorded;
- nearest-host assignment works on the initial frame.
