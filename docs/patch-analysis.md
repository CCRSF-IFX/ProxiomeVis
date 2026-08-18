# Patch Analysis

Patch Analysis displays precomputed results stored under
`uns["proxiome"]["patch"]` in the loaded H5AD. It does not call a patch
detection algorithm during an interactive session.

- **Markers** compares receiver and target marker frequencies and reports the
  stored patch-detection run plan.
- **Raji Signal** displays stored joint-proximity signal and the associated
  marker-abundance table.
- **Patch Burden** displays the stored per-cell patch-burden table.

The patch mapping can contain `run_plan`, `marker_unmixing`,
`raji_marker_abundance`, `raji_marker_proximity`, and `patch_burden` data
frames. Missing tables produce an unavailable-data message rather than an
inferred result.
