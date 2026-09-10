# Multipath Safety-Bound Gate

Date: 2026-09-10

## Question

Can a pairwise-certified grouped `c2` codebook fail only through collective
three/four-path dependence, and can a Gershgorin/coherence bound provide a
useful scalable safety certificate?

## Model

For a normalized finite-alphabet error event `delta`, define

```text
Phi_q(delta) = [H_1(q)delta, ..., H_P(q)delta]
G_q(delta) = Phi_q(delta)^H Phi_q(delta).
```

The exact conditioning metric is `lambda_min(G_q)`. For normalized path
columns, two sufficient lower bounds are

```text
lambda_min(G_q) >= min_i(1-sum_{j~=i}|G_q(i,j)|)
lambda_min(G_q) >= 1-(P-1)mu,
```

where `mu=max_{i~=j}|G_q(i,j)|`.

## Experiment

Script:

- `matlab_afdm/experiments/structural/run_multipath_safety_bound_gate.m`

Artifacts:

- `results/multipath_safety_bound_gate.csv`
- `results/multipath_safety_bound_by_event_class.csv`
- `results/multipath_safety_bound_gate.mat`
- `results/multipath_safety_bound_gate.png`

Configuration:

- `N=64`, `V=4`, integer Doppler, BPSK;
- path sets `(0,0),(2,2),(5,-3)`, `(0,0),(2,2),(8,0)`, and all four paths;
- baseline, rational GPS, local proposed, and pairwise-safe GPS patterns;
- 1652 events: 128 weight-one, 1000 sampled weight-two, 500 dense random,
  and 24 path-pair cycle events.

## Results

The Gershgorin bound is strongly correlated with the exact minimum eigenvalue:
approximately `0.9805~0.9913` across the twelve scheme/support combinations.
It gives a positive direct certificate for approximately `95.5%~99.9%` of
events.

However:

- all fourteen exact rank-loss occurrences are rational-GPS pair-cycle
  events;
- no collective-only rank loss is observed: every exact P-path rank loss is
  already detected as a pairwise dependence;
- the pairwise-safe GPS has no sampled rank-loss event for either P=3 or P=4;
- the minimum Gershgorin bound is negative for every scheme/support
  combination, so it cannot certify the complete worst-case event set without
  exact follow-up checks;
- the safe-but-uncertified fraction is about `0.1%~4.5%`;
- in the current MATLAB implementation and for tiny `P=3,4` Gram matrices,
  bound evaluation takes about `1.5~3.7` times as long as direct
  eigendecomposition because both require the Gram entries and the bound uses
  interpreted loops.

Representative exact minimum eigenvalues:

| Path set | Baseline | GPS rational | Local proposed | Pairwise-safe GPS |
|---|---:|---:|---:|---:|
| P3-set1 | 1.43e-4 | 0 | 1.77e-3 | 4.06e-2 |
| P3-set2 | 1.64e-1 | 5.55e-17 | 1.13e-1 | 4.27e-2 |
| P4-set3 | 1.43e-4 | 0 | 1.77e-3 | 3.24e-2 |

The pairwise-safe codebook improves the weakest tested P3-set1/P4-set3
conditioning substantially, but it is not uniformly best: baseline and local
proposed have larger margins in P3-set2.

## Decision

The initial random/event-class gate was NARROW, but subsequent directed search
changes the overall multipath decision to PURSUE as a thesis-scale candidate.

The bound is useful as an interpretable multipath diagnostic and as a hybrid
pre-filter that can certify most ordinary events. The current evidence does
not establish a new multipath-only failure mechanism or a computational
advantage strong enough for a separate core chapter.

To reopen this direction, first search specifically for a finite-alphabet
event that is pairwise safe but collectively rank deficient. If such events
remain absent under broader support/random searches, multipath analysis should
strengthen the first chapter rather than become the second chapter.

That reopening condition has now been met. An exhaustive valid `N=9` search
finds 564 strict collective-only events, and an `N=64` four-path construction
finds a pairwise-safe near-degenerate event with
`lambda_min=3.45e-10`. See
`research_notes/collective_multipath_rank_search.md`. The Gershgorin result
should now be treated as the first candidate certificate, not as the complete
solution.
