# Explicit c2 Side-Information Baseline

Date: 2026-09-10

## Purpose

Provide a minimal end-to-end interface for transferring the selected grouped
`c2` mode. This is a system-completeness baseline, not a thesis contribution.

## Implementation

- `afdm.tx.encode_mode_index` maps a one-based mode index to
  `ceil(log2(Q))` bits, applies configurable repetition, and BPSK modulates
  the coded bits.
- `afdm.rx.decode_mode_index` combines repeated BPSK observations and
  reconstructs the mode index.
- `run_explicit_c2_side_info_baseline` selects the minimum-PAPR mode from a
  four-mode grouped `c2` codebook, transmits the explicit mode header through
  an already-equalized AWGN control-channel abstraction, and uses the decoded
  mode for matched DAFT/MMSE data reception over a fractional doubly selective
  channel with perfect physical CSI.

The control-channel abstraction deliberately excludes synchronization and
channel estimation. It isolates mode-index errors and their propagation into
the AFDM payload.

## Representative run

- `N=64`, 16-QAM, four modes and two information bits per AFDM block;
- 1000 frames per data-SNR point;
- data SNR `[10,15,20]` dB;
- side-information SNR is 10 dB below data SNR;
- uncoded BPSK and threefold repetition.

Artifacts:

- `results/explicit_c2_side_info_baseline.csv`
- `results/explicit_c2_side_info_baseline.mat`
- `results/explicit_c2_side_info_baseline.png`

Results:

| Data/SI SNR | Method | Mode error | Data BER | Header symbols | Payload-relative overhead |
|---:|---|---:|---:|---:|---:|
| 10/0 dB | oracle | 0 | 0.1282 | 0 | 0 |
| 10/0 dB | BPSK-R1 | 0.167 | 0.1582 | 2 | 0.781% |
| 10/0 dB | BPSK-R3 | 0.016 | 0.1310 | 6 | 2.344% |
| 15/5 dB | oracle | 0 | 0.0488 | 0 | 0 |
| 15/5 dB | BPSK-R1 | 0.017 | 0.0521 | 2 | 0.781% |
| 15/5 dB | BPSK-R3 | 0 observed | 0.0488 | 6 | 2.344% |
| 20/10 dB | all | 0 observed | 0.00990 | 0/2/6 | 0/0.781/2.344% |

No-noise codec checks pass for all four modes and repetitions `1`, `3`, and
`5`.

## Interpretation

An incorrect mode causes a block-wide transform mismatch and increases data
BER, but a four-mode codebook needs only two information bits. Simple
repetition already makes the explicit header reliable at modest SNR for small
overhead. Therefore explicit side information is the preferred engineering
baseline; blind or joint mode detection must demonstrate a compelling benefit
before receiving further research effort.
