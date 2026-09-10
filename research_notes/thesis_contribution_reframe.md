# Thesis Contribution Reframe

Date: 2026-09-10

## System-completeness modules

The following are required for a credible end-to-end study but do not count as
independent contributions:

- explicit transfer and recovery of the selected `c2` mode;
- matched DAFT demodulation and MMSE reception;
- integer/fractional Doppler robustness checks;
- standard channel profiles and parameter justification;
- PAPR, BER/BLER, overhead and runtime reporting.

## Core work 1: exact two-path mechanism

Working title:

> Selection-aware diversity degradation and reliability-constrained
> pre-chirp codebook design for low-PAPR AFDM.

Contribution chain:

1. show that data-dependent `c2` selection invalidates analysis based on one
   fixed linear AFDM transform;
2. derive the two-path relative monomial/cycle condition;
3. connect exact alignment to rank loss and near alignment to finite-SNR
   coding-gain degradation;
4. construct and evaluate an offline safe codebook;
5. validate PAPR-reliability-complexity tradeoffs, including relative
   fractional Doppler.

## Core work 2 candidate: scalable multipath safety certification

The current closed-form mechanism is strongest for two paths. A general
`P`-path error event requires

```text
Phi_q(delta_x) = [H_1(q)delta_x, ..., H_P(q)delta_x]
```

to retain full column rank and acceptable conditioning over channel supports
and finite-alphabet error events. Pairwise safety does not automatically imply
full `P`-column safety, while exhaustive certification scales poorly.

Potential method:

- normalize path-response columns `u_p`;
- define worst-case mutual coherence
  `mu=max_{i~=j}|u_i^H u_j|`;
- use a Gershgorin-style sufficient bound
  `lambda_min(Phi^H Phi) >= 1-(P-1)mu` for normalized columns;
- construct a path-support conflict graph or robust codebook using this bound;
- compare analytical certification, exact numerical singular values and
  exhaustive small-system results;
- evaluate on three/four-path controlled channels and standard EVA/TDL
  profiles.

This direction matches a real neglected-structure and scalability gap rather
than inventing another PAPR optimizer. It becomes a core chapter only if the
bound is informative and enables a codebook that retains PAPR gain.

## Smallest gate for core work 2

1. Generate `P=3` and `P=4` path-response matrices for baseline, rational GPS,
   local proposed and integer-safe GPS patterns.
2. Compare the exact `lambda_min` with the coherence/Gershgorin lower bound.
3. Measure how often the bound is positive, how well it ranks codebooks, and
   whether it flags known dangerous patterns.
4. Compare certification cost with direct SVD over the same support/event set.

Decision:

- pursue if the bound is positive on a useful fraction of safe cases, detects
  dangerous cases, and reduces certification cost materially;
- narrow to a first-chapter analysis subsection if the bound is useful but too
  conservative for design;
- discard if it is almost always non-positive or does not preserve codebook
  ordering.

Initial random/event-class tests did not expose a collective-only loss.
Subsequent exhaustive `N=9` and directed `N=64` searches did: `N=9` contains
strict finite-alphabet three-path rank loss with every pair full rank, while
`N=64` contains a four-path pairwise-safe event with
`lambda_min=3.45e-10`. This establishes a real limitation of pairwise
certification and promotes scalable full-multipath safety design to the main
second-chapter candidate. See
`research_notes/collective_multipath_rank_search.md`.

## Lower-priority alternatives

- large-safe-codebook online branch-and-bound selection: only if the safe
  codebook must become large enough for exhaustive evaluation to be a real
  bottleneck;
- decoder/CRC-aided blind mode recovery: direct receiver coupling but narrow
  novelty because coded SLM precedents already exist;
- generic channel estimation, synchronization, MIMO and finite precision:
  currently outside the thesis main line.
