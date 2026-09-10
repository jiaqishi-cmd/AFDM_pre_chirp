# Dynamic c2 Common-Pilot Validation

## Question

Can a receiver use one mode-independent AFDM training block to estimate the physical channel, reconstruct the effective channel for every candidate grouped `c2` mode, and then identify the data-dependent mode under fractional Doppler?

This is a direction gate, not a complete receiver proposal.

## Controlled assumptions

- SISO AFDM, `N=32`, 16-QAM, three paths with delays `[0,1,2]`.
- Integer Doppler centers `[-2,0,2]` plus fractional offsets of magnitude `0`, `0.2`, or `0.4`.
- Pilot and data blocks share the same channel gains and lie in one coherence interval.
- Coarse delay/Doppler support is known; the fixed-`c2` training block estimates only complex path gains by LS.
- The data transmitter selects the lowest-PAPR mode from four representative four-group vector-`c2` patterns using local offsets `+-c2/16`.
- The receiver reconstructs a candidate-specific effective channel and uses an MMSE/slicer-residual metric for mode selection.

## Representative result

Run: 100 frames per point, SNR `[15,20,25]` dB, fractional magnitudes `[0,0.2,0.4]`.

Artifacts:

- `results/common_pilot_c2_validation_20260909_095440.csv`
- `results/common_pilot_c2_validation_20260909_095440.mat`
- `results/common_pilot_c2_validation_20260909_095440.png`

At 20 dB:

| Fractional magnitude | Gain NMSE | Mode accuracy | BER, estimated CSI and known mode | BER, estimated CSI and detected mode |
|---:|---:|---:|---:|---:|
| 0.0 | 8.61e-4 | 0.99 | 8.59e-3 | 9.77e-3 |
| 0.2 | 9.47e-4 | 0.90 | 1.56e-2 | 2.78e-2 |
| 0.4 | 9.08e-4 | 0.96 | 1.43e-2 | 1.91e-2 |

At 25 dB, mode accuracy is `0.96~1.00` and gain NMSE is approximately `2.6e-4~2.8e-4`.

## Interpretation

1. A common reference mode can provide a stable physical-gain estimate without knowing the data-selected `c2` mode.
2. Reconstructing all mode-specific effective channels from that common estimate is numerically consistent: effective-channel NMSE tracks path-gain NMSE.
3. Channel estimation is not the main loss in this controlled setting. Mode errors dominate the extra BER when the mode is unknown.
4. Fractional Doppler does not destroy the approach, but it changes mode-detection reliability even when the simple matrix-separation statistic remains nearly constant.
5. Therefore, Frobenius distance between candidate effective channels is not yet a sufficient codebook design metric. A receiver-aware metric should include the symbol alphabet, equalizer, channel uncertainty, and fractional leakage.

## Decision gate

The direction passes the first feasibility gate but is not yet a thesis contribution.

Continue only if the next compact study can establish both of the following:

- a mode-distance or pairwise-error metric that predicts the observed confusion behavior;
- acceptable mode accuracy after replacing known support with coarse delay/Doppler estimates and reducing the training overhead.

If either fails, fold fractional Doppler robustness into the transmitter safe-codebook chapter instead of building a separate receiver chapter.

## Follow-up 1: local-spacing sweep

The representative four-mode grouped codebook was swept over
`delta/c2 = [1/64, 1/32, 1/16, 1/8, 1/4]` at 20 dB and fractional
Doppler magnitude `0.4`, using 300 common random frames per point.

Artifacts:

- `results/common_pilot_c2_spacing_20260909_100556.csv`
- `results/common_pilot_c2_spacing_20260909_100556.mat`
- `results/common_pilot_c2_spacing_20260909_100556.png`

Key observations:

- Mean PAPR gain increases from `0.43 dB` at `1/64` to `1.01 dB` at
  `1/16`, then saturates near `1 dB`.
- Mode accuracy remains between `0.917` and `0.953`; larger geometric
  separation does not produce a monotonic accuracy improvement.
- Known-mode BER remains near `1.6e-2`, while detected-mode BER increases
  from `1.92e-2` to `2.75e-2` across the sweep.
- A larger spacing can make a wrong mode easier to distinguish in theory,
  but also makes the occasional wrong-mode decision more destructive.

This establishes a three-way tradeoff among PAPR gain, mode-error
probability, and mode-error consequence. Maximizing codebook distance alone
is not an appropriate objective.

## Follow-up 2: coarse Doppler-support mismatch

The true channel and receiver model were separated. The received pilot and
data use the true fractional Doppler, while the LS dictionary and reconstructed
effective channels use biased Doppler support.

Artifacts:

- `results/common_pilot_support_error_20260909_100942.csv`
- `results/common_pilot_support_error_20260909_100942.mat`
- `results/common_pilot_support_error_20260909_100942.png`

At 20 dB and fractional Doppler magnitude `0.4`:

| Support error (bin) | Effective-channel NMSE | Mode accuracy | Known-mode BER | Detected-mode BER |
|---:|---:|---:|---:|---:|
| 0 | 8.69e-4 | 0.947 | 1.59e-2 | 2.38e-2 |
| 0.02 | 1.69e-3 | 0.923 | 1.86e-2 | 3.01e-2 |
| 0.05 | 6.00e-3 | 0.893 | 2.93e-2 | 4.52e-2 |
| 0.10 | 2.13e-2 | 0.703 | 6.52e-2 | 1.05e-1 |
| 0.20 | 7.98e-2 | 0.460 | 1.53e-1 | 2.17e-1 |

The receiver is therefore sensitive to off-grid support mismatch. A viable
second thesis chapter cannot stop at LS gain estimation plus c2 detection. It
would need mode-aware continuous-Doppler refinement, joint mode/channel
updates, or another mechanism that explicitly controls this mismatch.

## Updated decision

The receiver direction remains viable, but its defensible focus is now
narrower:

> Dynamic grouped-pre-chirp AFDM reception with mode-aware off-grid Doppler
> refinement under imperfect channel support.

The next work should be analytical and algorithmic rather than a larger BER
run: formulate the joint residual in mode and continuous Doppler, then test a
small local-refinement step against the support-error curves above.
