# AFDM PAPR Prior-Art Map

Last updated: 2026-09-02

This note tracks papers and patents related to AFDM PAPR reduction, pre-chirp
parameter selection, side information, receiver recovery, and OFDM-style method
transfer. It is intended to protect novelty and guide thesis experiments.

## Current Project Position

Current working idea:

- Use the flexibility of AFDM pre-chirp parameter `c2`.
- Construct a small candidate set or group-wise candidate pattern.
- Select the candidate with lower PAPR.
- Evaluate PAPR, BER, receiver-side `c2` mismatch, semi-blind `c2` detection,
  and search complexity.

Main risk:

- Recent AFDM PAPR papers already use `c2` selection, SLM-like candidate
  generation, intrinsic/blind side information, affine-domain shifting, and
  low-complexity search.

Therefore, novelty should not be claimed as simply:

- "use `c2` to reduce PAPR"
- "select the minimum-PAPR AFDM candidate"
- "group subcarriers and vary pre-chirp parameters"

More defensible novelty should focus on the specific candidate construction,
receiver recovery mechanism, complexity reduction, robustness analysis, or a
clear application constraint.

## AFDM PAPR Papers Most Relevant To Novelty

### Yuan et al., PAPR Reduction with Pre-chirp Selection for AFDM

Source: https://arxiv.org/abs/2406.14064

Why it matters:

- This is the closest paper to the current idea.
- It proposes grouped pre-chirp selection (GPS).
- It varies the pre-chirp parameter across subcarrier groups and selects the
  lowest-PAPR candidate.
- It evaluates PAPR, complexity, spectral efficiency, and BER.

Novelty pressure:

- Any thesis or patent claim based only on group-wise `c2` selection will be
  weak against this paper.

Possible differentiation:

- Keep `c2` perturbations small around baseline rather than using GPS-style
  values.
- Emphasize structural-risk/BER robustness compared with GPS.
- Emphasize receiver-side candidate detection or intrinsic side information.
- Emphasize partial waveform reuse and topK complexity reduction.

### Nguyen and Bedeer, PAPR Reduction for AFDM by Affine-Domain Circular Shift without SI

Source: https://arxiv.org/html/2606.22464v1

Why it matters:

- Proposes affine-domain circular shift candidate generation.
- Selects the lowest-PAPR shifted candidate.
- Derives an ML-based receiver to detect the applied shift without side
  information using pilot/guard structure.
- Reports about 2.5 to 4 dB PAPR reduction.

Novelty pressure:

- Receiver-side blind/semi-blind candidate recovery is already a current
  research direction.

Possible differentiation:

- Compare your semi-blind `c2` candidate detector against this "no SI" framing.
- Avoid claiming blind recovery broadly unless the receiver evidence is strong.
- Focus on a simpler detector, smaller codebook, or compatibility with your
  pre-chirp perturbation pattern.

### Choi, Low-PAPR AFDM With Intrinsic SI and Low Complexity

Sources:

- https://www.semanticscholar.org/paper/Low-PAPR-AFDM-With-Intrinsic-SI-and-Low-Complexity-Choi/ee7c6a15cb9d0837d061ddc31fd9d2eea1584e3c
- https://www.scilit.com/ (search result text)

Why it matters:

- Proposes per-slot scalar chirp-offset selection.
- Claims low PAPR, intrinsic blind side information, and low receiver
  complexity.
- Search snippets describe a well-separated codebook where wrong candidates
  become slicer-inconsistent.

Novelty pressure:

- This overlaps with your current receiver-side `c2` detection experiment.

Possible differentiation:

- Deep-read this paper before finalizing the second thesis point.
- Identify whether your detector uses a different metric, candidate structure,
  channel-estimation assumption, or complexity path.

Evidence status:

- [INSUFFICIENT EVIDENCE] Full text not yet reviewed from an official accessible
  source in this note.

### Ali, Arous, and Arslan, Spreading the Wave: Low-Complexity PAPR Reduction for AFDM and OCDM

Source: https://arxiv.org/abs/2505.01778

Why it matters:

- Uses premodulation spreading to reduce PAPR in AFDM/OCDM.
- Compares WHT, DCT, Zadoff-Chu transform, and interleaved DFT.
- Positions itself as low-complexity and no-side-information compared with PTS
  and SLM.

Novelty pressure:

- Low-complexity and no-SI PAPR reduction are already explicit AFDM themes.

Possible differentiation:

- Use this as a baseline family or discussion point.
- Consider adding DFT/DCT/ZC spreading as baselines if implementation time
  allows.

### Reddy and Bitra, PAPR in AFDM: Upper Bound and Reduction With Normalized mu-Law Companding

Sources:

- https://ieeexplore.ieee.org/document/11002473/
- Cited and summarized in https://arxiv.org/html/2606.22464v1

Why it matters:

- Derives AFDM PAPR bounds and applies normalized mu-law companding.
- Reported as a low-complexity nonlinear waveform-domain PAPR method.

Novelty pressure:

- PAPR reduction with BER tradeoff is not new by itself.

Possible differentiation:

- Treat companding as a non-c2 baseline.
- Your method should report whether it avoids nonlinear distortion and what
  side-information/complexity cost it pays instead.

### Reddy and Bitra, PAPR Reduction in AFDM Using Zadoff-Chu Sequence Matrix and Companding Transform

Source: https://ieeexplore.ieee.org/document/11338808/

Why it matters:

- Uses Zadoff-Chu matrix precoding and companding transforms.
- Useful as a representative "precoding/sequence + companding" AFDM baseline.

Possible use:

- Include as literature comparison.
- Implement only if time allows.

### Karthiga and Deepa, Optimal Chirp Selection using Unimodular Quadratic Program

Source: https://www.sciencedirect.com/science/article/abs/pii/S187449072500134X

Why it matters:

- Uses optimization-style chirp selection.
- Relevant to search/optimization framing.

Evidence status:

- [INSUFFICIENT EVIDENCE] Abstract-level evidence only so far.

### Li et al., Chirp Parameter Selected Mapping for Low PAPR AFDM Transmission

Source: https://www.eurecom.edu/en/publication/8431

Why it matters:

- Explicitly frames AFDM PAPR reduction as chirp parameter selected mapping
  (CSM).
- Very close to "generate multiple c2 candidates, choose low PAPR".

Evidence status:

- [INSUFFICIENT EVIDENCE] Need full paper or official preprint.

### Cui et al., Adaptive c2-Perturbed AFDM Waveform Design for ISAC

Source: https://arxiv.org/abs/2606.04698

Why it matters:

- Directly optimizes subcarrier-wise `c2` perturbations.
- Jointly considers PAPR and autocorrelation sidelobes for ISAC.
- Uses closed-form gradients and spectral projected-gradient optimization.

Novelty pressure:

- "c2 perturbation" alone is no longer enough.

Possible differentiation:

- Your thesis can stay communication-centric instead of ISAC-centric.
- If expanding, add a small sensing/ACF metric only as optional future work,
  not a rushed main line.

### Gourar et al., Low-Complexity Sensing-Aware PAPR Reduction for AFDM-based ISAC Systems

Source: https://arxiv.org/html/2607.01064v1

Why it matters:

- Uses chirp-subcarrier reservation.
- Jointly reduces PAPR and improves ranging/sidelobe behavior.
- Uses gradient-based PAPR minimization plus randomized local search.

Possible use:

- Good source for "where AFDM PAPR is going" rather than an immediate thesis
  baseline.

## AFDM-Related Receiver And Channel Estimation Papers

### Zheng et al., Channel Estimation for AFDM With Superimposed Pilots

Sources:

- https://arxiv.org/abs/2404.10232
- https://pure.bit.edu.cn/en/publications/channel-estimation-for-afdm-with-superimposed-pilots/

Why it matters:

- Superimposed pilots improve spectral efficiency over guard-heavy embedded
  pilots.
- Iterative channel estimator and signal detector are relevant if your
  receiver-side `c2` detection assumes channel knowledge.

Use:

- Cite when discussing receiver assumptions and pilot overhead.

### Li et al., Matched Filtering-Based Channel Estimation for AFDM

Source: https://arxiv.org/abs/2507.09268

Why it matters:

- Low-complexity channel matrix construction and matched-filter channel
  estimation for doubly selective channels.
- Uses generalized Fibonacci search to reduce redundant computation.

Use:

- Relevant to your complexity story and receiver-side feasibility.

### Wu et al., BEM-Assisted Low-Complexity Channel Estimation for AFDM

Source: https://arxiv.org/abs/2504.18901

Why it matters:

- Uses basis expansion model to handle fractional Doppler and reduce exhaustive
  search.
- Connects estimation error to BER.

Use:

- Useful if your method is tested under estimated CSI instead of perfect CSI.

## Patent Landscape

### WO2022242870A1 / US20240094336A1, AFDM waveforms for doubly dispersive channels

Source: https://patents.google.com/patent/WO2022242870A1/

Why it matters:

- Broad AFDM waveform patent family.
- Covers generation/reception of multi-chirp AFDM waveforms, tunable
  coefficients, CSI/channel-parameter adaptation, pilot/guard designs, and
  signaling indications of parameters.

Novelty pressure:

- Avoid claiming the broad AFDM waveform or generic dynamic `c1/c2` adaptation.

### CN119094301A, 一种多维传输的峰均比抑制方法

Source: https://patents.google.com/patent/CN119094301A/zh

Why it matters:

- Chinese patent application for SLM-AFDM PAPR reduction.
- Builds multiple phase sequences, transforms them, calculates PAPR, and
  transmits the minimum-PAPR branch plus corresponding phase information.
- Cites Yuan's GPS paper.

Novelty pressure:

- Random phase-sequence SLM-style AFDM PAPR reduction is crowded.

Possible differentiation:

- Your method should avoid looking like generic SLM-AFDM.
- Emphasize `c2` perturbation structure, receiver-side candidate detection, and
  partial waveform reuse rather than arbitrary phase sequence multiplication.

### WO2024056176A1 / CN119856470A, spectral shaping signals in DAFT domain

Source: https://patents.google.com/patent/CN119856470A/zh

Why it matters:

- DAFT/AFDM-domain spectral shaping, pilot generation, and target spectral
  occupancy control.
- Not a direct PAPR candidate-selection patent, but it shows industrial
  patenting around DAFT-domain signal shaping and pilot design.

Use:

- Relevant to avoid overbroad claims about DAFT-domain shaping.

### OFDM PAPR patents

Representative source:

- https://patents.google.com/patent/CN107302516B/zh

Why it matters:

- Tone reservation, SLM, PTS, phase rotation, reserved carriers, and MIMO-OFDM
  PAPR suppression are mature patent areas.

Use:

- When transferring OFDM ideas to AFDM, novelty must come from AFDM-specific
  chirp/DAFT structure, receiver compatibility, or high-mobility channel
  behavior.

## OFDM Methods Worth Considering For AFDM Transfer

Likely useful:

- SLM-style selected mapping: already transferred to AFDM; high collision risk.
- PTS-style partitioning: possible but complexity-heavy; needs AFDM-specific
  partition/search reduction.
- DFT/DCT/WHT/ZC spreading: already explored for AFDM/OCDM; useful baseline.
- Tone/chirp-subcarrier reservation: active in AFDM-ISAC; useful if framed
  around low-complexity or sensing-aware metrics.
- Companding/clipping: easy baseline, but BER/OOB distortion must be reported.
- Active constellation extension: relevant but can complicate BER and receiver
  assumptions.
- Codebook design with intrinsic/blind SI: very relevant to your second point.

Less attractive under thesis time constraints:

- Full deep-learning PAPR reducers.
- Full ISAC waveform optimization.
- Heavy nonconvex optimization without a simple baseline implementation.
- Hardware PA modeling beyond a simple SSPA/Rapp model.

## Recommended Near-Term Search And Experiment Priorities

1. Deep-read Yuan GPS and Choi intrinsic-SI papers before finalizing claims.
2. Implement or approximate 2-3 baseline families:
   - GPS
   - DFT/DCT/ZC spreading
   - companding or SLM-AFDM
3. Strengthen your unique claim around:
   - small structured `c2` perturbation near baseline
   - receiver-side candidate detection
   - low-complexity partial waveform reuse/topK
   - robustness against GPS-style structural risk
4. Avoid rushing into ISAC unless advisor explicitly wants it.

