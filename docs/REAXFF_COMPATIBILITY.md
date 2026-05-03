# ReaxFF Compatibility

Compatibility outputs are preserved for downstream workflows:

- `fort.7` connection table
- `fort.8` connection table
- `molfra.out` per-frame molecule list / molecular fragmentation summary
- `xmolout`
- bond-order summaries when available
- charge summaries when available

Unless configured otherwise, these files are centroid-projected classical views
with cell-aware minimum-image topology analysis driven by the explicit cell in
the shared system state.  `fort.7` and `fort.8` are connection tables, not
generic snapshots.  `molfra.out` is a molecule-per-frame analysis with formula
counts and a component list.

PIMD-native bead distributions and estimators are written separately.

Current status: centroid-projected compatibility snapshots are now written by
the mixed BAOAB path at the end of a run.  They are intentionally separate from
the PIMD-native bead outputs and remain classical views of the mixed state.

The phase6-style temperature partition regression is a diagnostic gate only:
it measures carbon and electron temperature partitions and should not be read
as a claim that the hot-electron problem is physically solved in general.
