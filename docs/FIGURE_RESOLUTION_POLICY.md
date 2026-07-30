# Figure-resolution audit

The source raster figures in `05_figures` are audited with
`scripts/17_audit_figure_resolution.R`.

The default conservative gate is at least 600 effective dpi at a final width
of 180 mm. A per-file width can be supplied through the optional local file:

`06_logs_and_audit/figure_width_targets.tsv`

with columns `filename` and `target_width_mm`.

The script writes a detailed TSV and a human-readable report to
`06_logs_and_audit`. It does not modify or resave any figure. Word-embedded
copies can be downsampled by Word and are not a substitute for auditing the
separate source TIFF/PNG files intended for submission.
