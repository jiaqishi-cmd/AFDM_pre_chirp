# Curated Research Results

This directory contains small CSV and PNG artifacts selected for cross-machine
review. Large MAT files, logs, smoke tests, and intermediate outputs remain
ignored and can be reproduced from the committed MATLAB scripts.

Current groups:

- `explicit_c2_side_info_*`: explicit mode-index transmission baseline;
- `fractional_*` and `relative_fractional_*`: fractional-Doppler model and
  relative-Doppler closure checks;
- `mimo_c2_papr_*`: shallow MIMO separability gate;
- `multipath_safety_bound_*`: P=3/4 eigenvalue-bound diagnostics;
- `collective_only_rank_search_*`: exhaustive small-system collective-rank
  search;
- `n64_collective_rank_construction_*`: targeted N=64 four-path search.

Treat zero-error Monte Carlo entries as finite-sample observations, not proofs
of zero error probability. Read the corresponding files in `research_notes/`
before interpreting an artifact.
