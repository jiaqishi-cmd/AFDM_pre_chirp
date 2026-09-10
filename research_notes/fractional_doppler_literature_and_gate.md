# Fractional-Doppler AFDM Literature Map and Direction Gate

Date: 2026-09-09

## Purpose

Determine whether fractional Doppler creates a genuinely new failure mode for
selection-aware low-PAPR `c2` codebooks. Do not assume that fractional-Doppler
channel estimation or detection is the thesis direction.

## Core mechanism

Write each normalized Doppler shift as

```text
nu_p = alpha_p + beta_p, beta_p in (-0.5, 0.5].
```

For integer Doppler, one path contributes one cyclic shifted diagonal in the
DAFT-domain effective channel. Fractional `beta_p` changes the path response
into a Dirichlet-like leakage band. This has four consequences:

1. the exact monomial/cycle model used in the current integer-Doppler analysis
   no longer applies directly;
2. the effective channel is less sparse and low-complexity detectors retain
   more neighbours;
3. embedded pilots suffer inter-Doppler, inter-delay, pilot-pilot and
   pilot-data interference;
4. a larger guard parameter can capture more leakage but increases overhead.

The `c2` matrix is diagonal and unit-modulus, so it mainly changes effective
channel phases rather than the magnitude envelope of fractional-Doppler
leakage. The thesis-relevant question is therefore not whether `c2` changes
the leakage width, but whether leakage changes the rank, singular values and
PEP ordering of data-dependent `c2` codebooks.

## Primary literature

### 1. AFDM foundation and fractional input-output model

- A. Bemani, N. Ksairi, M. Kountouris, "Affine Frequency Division Multiplexing
  for Next Generation Wireless Communications," IEEE TWC, 2023.
- https://arxiv.org/abs/2204.12798

Useful content:

- exact integer/fractional DAFT-domain path response;
- `nu_p = alpha_p + beta_p` decomposition;
- guard parameter in the `c1` setting;
- embedded-pilot estimation that first locates integer delay/Doppler and then
  searches/refines the fractional component;
- full-diversity proof is explicit for integer Doppler, while the fractional
  extension is stated with modifications and supported numerically.

### 2. Independent doubly-dispersive MATLAB reference

- H. S. Rou et al., "From OTFS to AFDM: A Comparative Study of Next-Generation
  Waveforms for ISAC in Doubly-Dispersive Channels," IEEE SPM, 2024.
- Paper: https://arxiv.org/abs/2401.07700
- Code: https://github.com/eric-hs-rou/doubly-dispersive-channel-simulation

Reusable files:

- `AFDMmod.m`, `AFDMdemod.m`;
- `cconvLTVChannel.m`, which accepts continuous Doppler frequencies;
- `sampleCode.m` and the effective-channel visualizer.

Use this package as an independent reference for channel-matrix and AFDM
transform conventions. Do not copy it into the production source tree before
the model comparison passes.

### 3. Fractional-Doppler channel estimation by diagonal reconstruction

- H. Yin et al., "Diagonally Reconstructed Channel Estimation for MIMO-AFDM
  With Inter-Doppler Interference in Doubly Selective Channels," IEEE TWC,
  2024.
- https://arxiv.org/abs/2206.12822

Useful content:

- explicitly separates IDoI, IDI, IPI and IPDI;
- shows that estimating every physical path can be fragile when paths share a
  delay and have different fractional Dopplers;
- reconstructs the effective channel directly from a pilot column;
- reports the guard-width versus overhead tradeoff and a non-monotonic pilot
  power effect.

Borrow only its leakage-band, retained-energy and truncation diagnostics in
the first stage. EPA-DR itself is a receiver baseline, not the current thesis
direction.

### 4. Fractional-Doppler channel estimation without known path count

- H. Jia et al., "Pilot-based PSC Fractional Channel Estimation for
  Self-Adaptive AFDM Communications," APCC 2025.
- https://www.ieice.org/publications/proceedings/summary.php?expandable=0&iconf=APCC&number=W2-1-4&session_num=W2-1&year=2025

It transforms the DAFT observation into a pseudo delay-Doppler domain and uses
pilot/checking-function correlation peaks. It is useful for understanding
off-grid visualization and peak-search refinement, but is currently outside
scope.

### 5. Fractional delay/Doppler local refinement

- "Joint Fractional Delay and Doppler Frequency Estimator Under Spectrum
  Wrapping Phenomenon for LEO-ICAN AFDM Signals," 2026.
- https://arxiv.org/abs/2602.04316

It combines peak-to-sidelobe-power-ratio detection with an early-late gate.
The early-late local interpolation idea is a possible lightweight reference if
the waveform study later requires continuous-Doppler refinement.

### 6. Detection under leakage

- L. Wu et al., "AFDM Signal Detection Based on Message Passing Scheme,"
  Digital Signal Processing, 2024.
- https://doi.org/10.1016/j.dsp.2024.104633

The ordinary MP detector loses its sparse-graph advantage under fractional
Doppler. The proposed MF-MP detector uses matched-filter diagonal enhancement.
Useful diagnostic: compare sparsity/energy concentration of `H_eff` and
`H_eff^H H_eff`.

- C.-W. Chen et al., "PCG-Based Expectation Propagation Detector With
  Message-Passing Initialization and Approximated Variance for Fractional
  Doppler AFDM," IEEE WCL, 2025.
- https://doi.org/10.1109/LWC.2025.3627115

This is a strong recent detection baseline, but implementing it now would move
the project away from the existing contribution.

### 7. Pilot-overhead alternatives

- GI-free AFDM channel estimation: https://arxiv.org/abs/2404.01088
- Superimposed pilots: https://arxiv.org/abs/2404.10232
- Single-pilot interference-position method:
  https://doi.org/10.1109/LWC.2025.3573104

These papers show that fractional-Doppler receiver and pilot design is already
crowded. They should be literature baselines only unless the thesis is later
deliberately changed to receiver research.

## What to test first

### Experiment F0: independent model agreement

Compare the current `afdm.channel.multipath` and
`afdm.rx.estimate_effective_channel` against Rou's
`cconvLTVChannel.m` for identical path gains, delays and fractional Dopplers.

Metrics:

- relative time-domain channel-matrix error;
- relative DAFT-domain effective-channel error;
- integer and fractional test cases;
- explicit sign/convention audit for Doppler and CPP.

Kill condition: do not interpret later results until both models agree after
documented convention conversion.

### Experiment F1: one-path leakage anatomy

Sweep `beta = [0, 0.1, 0.2, 0.3, 0.4, 0.5]` for one path.

Measure:

- effective-channel heat map;
- energy captured by top `K` entries per row;
- 90%, 95% and 99% leakage bandwidth;
- approximation error after band/top-K truncation;
- dependence on `N`, `c1` guard parameter and path delay.

Do not treat `c2` as a leakage-magnitude control variable without evidence;
unit-modulus pre-chirp phases normally preserve entry magnitudes.

### Experiment F2: exact fractional pairwise-diversity audit

Replace integer monomial path matrices in the existing pairwise audit with the
exact fractional path matrices. For baseline, GPS, proposed and safe-GPS
patterns, sweep one or two fractional components over `[-0.5, 0.5]`.

Measure:

- rank of `Phi(delta)`;
- minimum non-zero singular value;
- product distance / determinant surrogate;
- PEP bound;
- whether the integer dangerous point becomes regularized, remains poorly
  conditioned, or moves to another fractional value.

This is the decisive experiment for thesis relevance.

### Experiment F3: perfect-CSI BER confirmation

Only after F2 identifies representative safe and unsafe fractional cases, run
small perfect-CSI ML/MMSE BER comparisons. Exclude channel estimation so the
result isolates waveform/codebook geometry.

Baselines:

- conventional AFDM;
- Yuan GPS;
- current local-perturbation method;
- integer-safe codebook;
- a fractional-robust candidate if F2 supports one.

### Experiment F4: guard/performance tradeoff

Sweep the fractional guard parameter after a waveform effect is established.
Report BER/PEP improvement together with retained bandwidth, effective sparsity
and spectral/pilot overhead. Do not present a larger guard as a free gain.

## Thesis decision rule

PURSUE a second chapter only if all of the following are observed:

1. the integer-Doppler safety criterion fails or becomes materially inaccurate
   under fractional Doppler;
2. the failure has a reproducible mechanism in singular-value/PEP analysis;
3. a new leakage-aware robust constraint or codebook construction restores
   reliability over a fractional-Doppler interval;
4. the robustness cost in PAPR, complexity and overhead is measurable and
   acceptable.

If fractional Doppler only changes channel-estimation or detector performance,
while the relative safety ordering of the `c2` codebooks remains unchanged,
keep it as a robustness subsection of the first core chapter rather than a
second chapter.

## Current recommendation

Run F0-F2 only. Do not implement PSC, EPA-DR, EP, MF-MP, GI-free pilots or the
previous common-pilot receiver until the waveform-level gate is passed.

## Initial F0-F2 results

### F0 model cross-check

Script:

- `matlab_afdm/experiments/structural/run_fractional_model_crosscheck.m`

Artifacts:

- `results/fractional_model_crosscheck.csv`
- `results/fractional_model_crosscheck.mat`

The current channel implementation agrees with a separate explicit
sample-wise complex-exponential construction to approximately `4e-16` for
all tested fractional Dopplers. The downloaded reference package agrees at
integer Doppler, but differs off grid:

| beta | Current vs explicit | Current vs downloaded reference |
|---:|---:|---:|
| 0 | 3.67e-16 | 3.79e-16 |
| 0.10 | 3.78e-16 | 2.99e-1 |
| 0.25 | 3.64e-16 | 6.85e-1 |
| 0.40 | 3.90e-16 | 9.21e-1 |
| 0.50 | 3.61e-16 | 9.68e-1 |

The discrepancy is consistent with the downloaded `cconvLTVChannel.m` using
the non-integer matrix power `W^nu`, whose principal-branch convention is not
the same as the unwrapped sample-wise phase `exp(-j*2*pi*nu*n/N)`. Therefore,
the package is useful for integer-model cross-checks and visualization, but
its fractional path generator must not be treated as an oracle without a
documented correction.

### F2 c2 safety gate

Script:

- `matlab_afdm/experiments/structural/run_fractional_doppler_c2_safety_gate.m`

Configuration:

- `N=64`, `V=4`, `c1=7/(2N)`;
- BPSK candidate error set with all weight-one events, 1500 deterministic
  weight-two samples, and each scheme's integer cycle event;
- second-path fractional Doppler swept over `[-0.5,0.5]`;
- schemes: baseline, rational GPS pattern `[2 2 1 1]`, local proposed pattern,
  and phase-offset integer-safe GPS;
- three known integer-danger supports: `(delay,Doppler)=(2,2),(5,-3),(8,0)`.

Artifacts use tags `d2_a2`, `d5_a-3`, and `d8_a0`, for example:

- `results/fractional_doppler_c2_safety_gate_d2_a2.csv`
- `results/fractional_doppler_c2_safety_gate_d2_a2.png`

Minimum singular values at `beta=0` and over the sampled nonzero fractional
points are:

| Support | Scheme | beta=0 | Minimum off grid |
|---|---|---:|---:|
| (2,2) | baseline | 0.442 | 0.445 |
| (2,2) | GPS-rational | 0 | 0.0641 |
| (2,2) | local-proposed | 0.336 | 0.341 |
| (2,2) | integer-safe-GPS | 0.296 | 0.302 |
| (5,-3) | baseline | 0.0120 | 0.0652 |
| (5,-3) | GPS-rational | 0 | 0.0641 |
| (5,-3) | local-proposed | 0.0420 | 0.0766 |
| (5,-3) | integer-safe-GPS | 0.211 | 0.220 |
| (8,0) | baseline | 0.558 | 0.560 |
| (8,0) | GPS-rational | 0 | 0.0636 |
| (8,0) | local-proposed | 0.533 | 0.534 |
| (8,0) | integer-safe-GPS | 0.209 | 0.216 |

At `|beta|=0.5`, one DAFT-domain path retains only about 40.5% of its energy
in its largest entry, 85.6% in its largest three entries, and needs about
eight entries for 95% energy. Across all four `c2` profiles, the maximum
difference in path-matrix entry magnitudes is below `4.5e-16`, confirming that
`c2` changes phase geometry rather than the leakage-magnitude envelope.

### Interpretation

The first evidence does not support promoting fractional Doppler to a second
core chapter:

1. the exact rational-GPS cycle loss occurs at the integer point;
2. moving off grid regularizes that exact rank loss rather than creating a new
   one in the tested event set;
3. the local proposed and integer-safe GPS profiles remain non-singular over
   the tested fractional interval;
4. the integer-safe GPS margin does not collapse at the three known dangerous
   supports;
5. leakage strongly affects channel sparsity and receiver complexity, but its
   magnitude is effectively independent of `c2`.

There is still a near-integer vulnerability: at the first sampled nonzero
offset, rational GPS has minimum sigma near `0.064`, so the integer danger does
not disappear instantly. However, this currently strengthens the first
chapter's finite-SNR conditioning argument; it does not yet require a new
fractional-robust codebook.

The experiment is not exhaustive for all weight-three-or-higher fractional
error events or all joint fractional components. A broader support/event audit
could change the conclusion, but it should be run only if needed to establish
robustness of the first chapter, not on the assumption that a second chapter
already exists.

## Final relative-Doppler closure

Script:

- `matlab_afdm/experiments/structural/run_relative_fractional_doppler_closure.m`

Artifacts:

- `results/relative_fractional_doppler_closure.mat`
- `results/relative_fractional_doppler_closure.png`
- `results/relative_fractional_doppler_ber.csv`
- `results/relative_fractional_doppler_ber.png`

For the rational-GPS integer danger event at support `(2,2)`, both path
fractional components were swept over a two-dimensional grid. Points with the
same relative fractional Doppler `delta_beta=beta2-beta1` produce the same
minimum singular value to within `1.11e-15`. The complete diagonal
`beta1=beta2` remains rank deficient, including `(0.25,0.25)` and `(0.5,0.5)`.

Perfect-CSI fixed-channel MMSE simulations with 3000 frames per point confirm
the geometry:

- at SNR 28 and 40 dB, rational GPS retains a BER near `4e-3` for common
  fractional components `0`, `0.25`, and `0.5`;
- baseline, local proposed, and integer-safe GPS record no errors at these
  aligned points in the sampled frames;
- with relative fractional offset `0.05`, the rational-GPS BER falls to about
  `5.2e-6` at 28 dB and no sampled errors at 40 dB;
- with relative offset `0.25`, all schemes record no sampled errors.

Final statement for the first chapter:

> The GPS risk is governed by relative delay-Doppler cycle alignment rather
> than by whether each absolute Doppler is integer. A common fractional
> component cancels in the relative path operator and preserves the exact
> cycle; a nonzero relative fractional offset breaks the exact cycle, while a
> small offset leaves a finite-SNR conditioning neighbourhood.

Fractional Doppler is now closed as a robustness subsection of the first core
chapter. No additional pilot, channel-estimation, or fractional detector work
is justified by the waveform-level evidence.
