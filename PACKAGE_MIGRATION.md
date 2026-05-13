# MATLAB Package Migration Plan

Goal: gradually move stable code from global path-based functions to MATLAB
package namespaces. This reduces name collisions and lets future code expose
only the `matlab_afdm/` root instead of every subfolder.

## Long-Term Target

Recommended final structure:

```text
matlab_afdm/
  +afdm/
    +config/
    +tx/
    +rx/
    +channel/
    +chirp/
    +metrics/
    +search/
    +analysis/
    +delta/
  experiments/
```

New code should prefer calls like:

```matlab
metrics = afdm.metrics.c2_structural(c2Vec, cfg);
profile = afdm.chirp.build_profile('proposed_grouping', N, params);
```

## MATLAB Constraint

MATLAB packages still require the package parent directory to be visible.

Practical meaning:

- If current folder is `matlab_afdm/`, `afdm.*` calls work directly.
- If running from elsewhere, only `matlab_afdm/` needs to be on the path.
- We no longer need to expose every implementation subfolder once migration is complete.

## Current Migration Status

### Shared analysis and delta helpers

Done:

- Added `matlab_afdm/+afdm/+analysis/`
- Added:
  - `afdm.analysis.build_path_matrix`
  - `afdm.analysis.bemani_equation_error`
  - `afdm.analysis.phi_metrics`
  - `afdm.analysis.dmin_over_delta_set`
- Added `matlab_afdm/+afdm/+delta/`
- Added:
  - `afdm.delta.generate_set_extended`
  - `afdm.delta.generate_set_for_dmin`
  - `afdm.delta.generate_set_for_equation_search`
  - `afdm.delta.build_recursive_from_gps`
- Removed root-level analysis/delta helper wrappers after updating call sites.
- Moved formal reproduction and resume runners to
  `matlab_afdm/experiments/reproduction/`.

### Stage 1: metrics package

Done:

- Added `matlab_afdm/+afdm/+metrics/`
- Added:
  - `afdm.metrics.c2_structural`
  - `afdm.metrics.phase_structural`
  - `afdm.metrics.c2_candidate_correlation`
- Kept compatibility wrappers:
  - `calc_c2_structural_metrics`
  - `calc_phase_structural_metrics`
  - `calc_c2_candidate_correlation`
- Updated selected experiments to call `afdm.metrics.*`.

### Stage 2: chirp package

Done:

- Added `matlab_afdm/+afdm/+chirp/`
- Added:
  - `afdm.chirp.build_profile`
  - `afdm.chirp.baseline_profile`
  - `afdm.chirp.gps_profile`
  - `afdm.chirp.proposed_profile`
  - `afdm.chirp.build_gps_pattern`
  - `afdm.chirp.build_proposed_pattern`
  - `afdm.chirp.group_index`
  - `afdm.chirp.gps_candidate_set`
  - `afdm.chirp.get_param`
- Updated core pre-chirp entry points:
  - `afdm.chirp.apply_scheme`
  - `afdm.chirp.select_for_symbols`
  - `afdm.chirp.select_greedy_profile`
- Removed local compatibility wrappers after archiving legacy code:
  - `build_pre_chirp_profile`
  - `build_c2m_gps_pattern`
  - `build_c2m_proposed_pattern`
  - `pre_chirp/`

### Stage 3: search package

Partially done:

- Added `matlab_afdm/+afdm/+search/`
- Added:
  - `afdm.search.greedy_group_papr_selection`
  - `afdm.search.full_beam_search`
  - `afdm.search.reuse_beam_search`
  - `afdm.search.precompute_partial_waveforms`
  - `afdm.search.combine_partial_waveform`
  - `afdm.search.direct_full_waveform`
  - `afdm.search.full_waveform`
  - `afdm.search.ifft_oversampled`
- Kept compatibility wrapper:
  - `greedy_group_papr_selection`
- Updated:
  - `select_greedy_profile`
  - `experiments/complexity/run_partial_reuse_theory_and_timing.m`

Remaining search cleanup:

- Review whether the package-level `afdm.search.*` functions should expose a
  smaller public facade for experiments, or remain as implementation helpers.

### Stage 4: tx/rx/channel packages

Stage 4A tx package is done:

- Added `matlab_afdm/+afdm/+tx/`
- Added:
  - `afdm.tx.engine`
  - `afdm.tx.idaft_mod`
  - `afdm.tx.add_cpp`
  - `afdm.tx.compute_papr`
  - `afdm.tx.random_data`
- Removed local compatibility wrappers after archiving legacy code:
  - `afdm_tx_engine`
  - `idaft_mod`
  - `add_cpp`
  - `compute_papr`
  - `random_data_generator`
- Updated:
  - `simulate_frame`
  - `afdm.search.*` PAPR/IDAFT calls

Stage 4B rx/channel package is done:

- Added `matlab_afdm/+afdm/+rx/`
- Added:
  - `afdm.rx.engine`
  - `afdm.rx.remove_cpp`
  - `afdm.rx.daft_demod`
  - `afdm.rx.estimate_effective_channel`
  - `afdm.rx.mmse_equalize`
  - `afdm.rx.equalize_symbols`
  - `afdm.rx.symbol_decision`
  - `afdm.rx.compute_bit_errors`
- Added `matlab_afdm/+afdm/+channel/`
- Added:
  - `afdm.channel.generate_profile`
  - `afdm.channel.multipath`
  - `afdm.channel.add_awgn`
- Removed local compatibility wrapper directories after archiving legacy code
  to remote branches:
  - `legacy/pre-package-migration`
  - `legacy/pre-rx-channel-package`
- Updated:
  - `afdm_config`
  - `configure_experiment`
  - `simulate_frame`
  - channel diagnostic helpers

Remaining package cleanup:

- Review older exploratory scripts for possible archival, but packageable shared
  helpers have been moved under `+afdm`.

### Stage 5: reduce setup_paths

`setup_paths` now exposes only `matlab_afdm/` and `experiments/`.

## Compatibility Policy

The current branch has moved past the local-wrapper compatibility stage for
tx/rx/channel/chirp. Older path-based code is preserved on the remote legacy
branches, while active code should call `afdm.*` package functions directly.
