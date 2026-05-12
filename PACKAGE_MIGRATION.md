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
- Kept compatibility wrappers:
  - `build_pre_chirp_profile`
  - `build_c2m_gps_pattern`
  - `build_c2m_proposed_pattern`
- Updated core pre-chirp entry points:
  - `apply_pre_chirp_scheme`
  - `select_greedy_profile`

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

- Update `run_partial_reuse_M_sweep.m`
- Update `run_partial_reuse_topK_sweep.m`
- Consider migrating `run_papr_search_complexity_study.m`

## Recommended Next Steps

### Stage 4: tx/rx/channel packages

Move the core link functions into:

```text
+afdm/+tx/
+afdm/+rx/
+afdm/+channel/
```

Candidates:

- `idaft_mod`, `compute_papr`, `afdm_tx_engine`
- `daft_demod`, `mmse_equalize`, `estimate_effective_channel`
- `multipath_channel`, `add_awgn`

This stage is riskier because many functions call each other. Do it in a
separate commit with smoke tests.

### Stage 5: reduce setup_paths

Once all stable functions are packaged, `setup_paths` can stop adding every
implementation subfolder. It should only expose `matlab_afdm/` and, if useful,
`experiments/`.

## Compatibility Policy

For each migrated module, keep the old function name as a wrapper for at least
one migration stage:

```matlab
function y = old_function(varargin)
    y = afdm.module.new_function(varargin{:});
end
```

After all scripts switch to package calls and smoke tests pass, wrappers can be
removed in a dedicated cleanup commit.
