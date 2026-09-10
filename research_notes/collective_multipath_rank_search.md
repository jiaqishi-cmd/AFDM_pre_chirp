# Collective Multipath Rank Search

Date: 2026-09-10

## Objective

Find a finite-alphabet event for which every two-path response is full rank,
but the complete P-path response matrix is rank deficient or severely ill
conditioned. Such an event proves that pairwise codebook certification is not
sufficient for general multipath reliability.

## Exhaustive small-system search

Script:

- `matlab_afdm/experiments/structural/run_collective_only_rank_search.m`

The valid exhaustive configuration is:

- `N=9`, `V=3`, `alpha_max=1`, `max_delay=2`;
- all `3^9-1=19682` nonzero BPSK difference vectors;
- nine integer delay-Doppler supports and all 84 three-support subsets;
- baseline, all eight rational-GPS patterns, all eight phase-offset GPS
  patterns, and all 27 local-perturbation patterns.

Results:

| Family | Pairwise-safe patterns | Patterns with collective-only loss | Collective-only events |
|---|---:|---:|---:|
| baseline | 1/1 | 0 | 0 |
| GPS-rational | 8/8 | 8 | 384 |
| GPS-offset | 8/8 | 8 | 180 |
| local-proposed | 27/27 | 0 | 0 |

A representative rational-GPS event uses paths `(0,-1)`, `(1,-1)`, and
`(2,-1)`, and a weight-two BPSK difference at positions 2 and 5. Direct SVD
gives

```text
P=3 singular values = [1.22474487139, 1.22474487139, 1.46e-15]
rank = 2
```

Every two-column submatrix has singular values

```text
[1.22474487139, 0.707106781187]
```

and pairwise correlation magnitude `0.5`. The three-column null residual is
`1.46e-15`. This is a strict collective-only finite-alphabet rank-loss event,
not a tolerance artefact.

A legal power-of-two check with `N=8`, `alpha_max=0`, and delays `0:7` finds no
collective-only event in all 6560 BPSK differences. Thus exact collective loss
depends on transform/cycle arithmetic and is not universal at every size.

Artifacts:

- `results/collective_only_rank_search_N9_V3_a1_l2.mat`
- `results/collective_only_rank_search_N9_V3_a1_l2_findings.csv`
- `results/collective_only_rank_search_N9_V3_a1_l2.png`
- `results/collective_only_rank_search_N8_V4_a0_l7.mat`

## N=64 targeted four-path construction

Script:

- `matlab_afdm/experiments/structural/run_n64_collective_rank_construction.m`

The four supports `(0,0)`, `(2,2)`, `(5,-3)`, and `(7,-1)` have effective
integer shifts near `0`, `16`, `32`, and `48`. A structured event set places
ternary BPSK differences at equal offsets in four 16-sample blocks. All 1280
events and 114 grouped patterns are evaluated.

No exact collective-only rank loss is found. However, the local-proposed
pattern `[1,3,1,1]` has a pairwise-independent weight-two event with

```text
event positions = [32,48]
block amplitudes = [0,-1,-1,0]
P=4 singular values = [1.41421356225, 1.00001313692,
                       0.99998686291, 1.85785287568e-05]
condition number = 7.61e4
lambda_min(Phi^H Phi) = 3.45e-10
```

Every two-path submatrix remains well conditioned:

```text
minimum pair singular value = 0.707106781187
maximum pair correlation = 0.5
```

Therefore pairwise screening misses a severe P=4 collective conditioning
event at the actual `N=64` system size. It is not exact rank loss, but it can
produce a near-indistinguishable finite-alphabet error event under a matching
four-path gain vector.

Family-level minimum pairwise-independent eigenvalues in this targeted set:

| Family | Best pattern | Worst pattern |
|---|---:|---:|
| baseline | 3.57e-5 | 3.57e-5 |
| GPS-rational | 8.93e-2 | 8.93e-2 |
| GPS-offset | 4.19e-2 | 1.61e-2 |
| local-proposed | 2.80e-4 | 3.45e-10 |

Artifacts:

- `results/n64_collective_rank_construction.mat`
- `results/n64_collective_rank_construction.csv`
- `results/n64_collective_rank_construction.png`

## Updated decision

PURSUE as a thesis-scale candidate, subject to one more occurrence/protection
gate.

The neglected structure is now concrete:

> Pairwise-safe grouped pre-chirp patterns can still be collectively rank
> deficient or severely ill conditioned when three or more path responses are
> considered jointly.

The next method should target the full Gram matrix rather than only path
pairs. Candidate approaches include a tighter multi-column eigenvalue bound,
incremental Cholesky/Schur-complement certification, and codebook construction
that maximizes a worst-case P-path conditioning margin.

Before promoting this to the final second chapter, verify that the N=64 near
event produces a PEP/BER penalty and determine how often PAPR selection chooses
patterns with poor collective margins. A purely adversarial event with zero
selection probability would remain a useful robustness result but not a full
chapter motivation.
