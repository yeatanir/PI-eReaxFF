# ReaxFF Compatibility

Compatibility outputs are preserved for downstream workflows:

- `fort.7`
- `fort.8`
- `molfra.out`
- `xmolout`
- bond-order summaries when available
- charge summaries when available

Unless configured otherwise, these files are centroid-projected classical views.
PIMD-native bead distributions and estimators are written separately.

Current status: centroid-projected compatibility snapshots are now written by
the mixed BAOAB path at the end of a run.  They are intentionally separate from
the PIMD-native bead outputs and remain classical views of the mixed state.
