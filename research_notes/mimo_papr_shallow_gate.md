# MIMO c2 PAPR Shallow Gate

Date: 2026-09-10

## Question

Does extending the SISO `c2` PAPR selector to multiple independent transmit
antennas create a nontrivial joint minimax-PAPR optimization problem?

## Separable model

For antenna `t`, let `f_t(q_t)` denote the PAPR obtained with mode `q_t`. If
every antenna may select its mode independently and there are no spatial-code,
precoder, shared-mode, or PA-coupling constraints, then

```text
min_{q1,...,qNt} max_t f_t(q_t) = max_t min_{qt} f_t(q_t).
```

Therefore the globally optimal Cartesian search is achieved by independently
selecting the minimum-PAPR mode on every antenna. The apparent `Q^Nt` search
space does not by itself create a new algorithmic problem.

## Experiment

Script:

- `matlab_afdm/experiments/papr/run_mimo_c2_papr_separability_gate.m`

Artifacts:

- `results/mimo_c2_papr_separability_gate.csv`
- `results/mimo_c2_papr_separability_gate.mat`
- `results/mimo_c2_papr_separability_gate.png`

Configuration:

- `N=64`, 16-QAM, 5000 frames;
- four representative grouped local-perturbation `c2` modes;
- transmit-antenna counts `2`, `4`, and `8`;
- compare independent per-antenna selection, exhaustive Cartesian minimax
  selection, and a constrained common-mode selector.

Results:

| Nt | Cartesian combinations | Independent mean worst PAPR | Common-mode mean PAPR | Mean common penalty | Independent selections all equal |
|---:|---:|---:|---:|---:|---:|
| 2 | 16 | 6.168 dB | 6.432 dB | 0.264 dB | 24.68% |
| 4 | 256 | 6.473 dB | 7.000 dB | 0.526 dB | 1.74% |
| 8 | 65536 | 6.748 dB | 7.516 dB | 0.768 dB | 0% sampled |

For `Nt=2` and `Nt=4`, exhaustive Cartesian minimax search and independent
selection agree on every frame with maximum numerical error zero. Exhaustive
enumeration is omitted for `Nt=8` because the analytical identity already
establishes the result.

The common-mode restriction reduces mode signalling from `Nt*log2(Q)` to
`log2(Q)` bits, but produces a mean PAPR penalty that grows with antenna count;
the 99th-percentile penalty is approximately `1.56`, `1.76`, and `2.00` dB for
`Nt=2`, `4`, and `8`.

## Decision

Plain MIMO extension is not a second thesis contribution. It becomes
nonseparable only after introducing at least one substantive coupling:

- one common `c2` mode or restricted cross-antenna mode tuples;
- STBC or another spatial codeword constraint;
- precoding/beamforming that mixes antenna signals;
- a spatial-diversity/product-distance constraint;
- shared or coupled PA constraints.

Each of these requires a new MIMO system model and substantial receiver or
codeword analysis. Given the current timeline, keep MIMO as a low-priority
candidate and do not build the full MIMO-AFDM chain unless literature review
reveals a sharply bounded, unoccupied coupling problem.
