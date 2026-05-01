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

Current status: ReaxFF compatibility output is still an explicit stub.  The
harmonic driver currently writes only PIMD-native outputs, which keeps the
scientific outputs separate from future centroid-projected compatibility files.
