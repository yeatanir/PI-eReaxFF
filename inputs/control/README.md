# Control files for `pair_style reaxff/metal <control> ereaxff yes ...`
The control parser accepts no comments. Two files:
- `legacy.ctrl` — the legacy thresholds (identical to `NULL`, verified). Use for ALL validation / reference comparisons against `ereac`.
- `smooth.ctrl` — nbrhood_cutoff 5.2, thb_cutoff 1e-9, thb_cutoff_sq 1e-18. Use for any run with MOBILE atoms. Removes the two energy jumps of `validation/phase3_engine/DISCONTINUITY.md`; changes the energy surface (`validation/phase3_engine/SMOOTH_ADOPTION.md`).
Every run script logs the control file it used (`CONTROL_USED` line / `control_used.txt`).
