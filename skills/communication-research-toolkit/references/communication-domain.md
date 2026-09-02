# Communication Domain Checks

Use these checks when reading, reviewing, or extending communication research.

System model:
- Identify modulation, frame structure, channel model, mobility model, synchronization assumptions, and receiver knowledge.
- Check whether the paper assumes perfect CSI, estimated CSI, partial CSI, or no CSI.
- Check whether the model includes Doppler, delay spread, carrier frequency offset, timing offset, phase noise, hardware impairment, or PAPR constraints.

AFDM and waveform research:
- Track DAFT/IDAFT definitions, chirp parameters, guard design, pilot design, and receiver structure.
- Separate waveform-domain gains from receiver, channel, or simulation-setting gains.
- Check whether PAPR, BER, spectral efficiency, complexity, and robustness are all considered when relevant.

Experimental fairness:
- Confirm that baselines use comparable power, bandwidth, coding, modulation, channel realization, receiver knowledge, and search budget.
- Treat gains as weak if the proposed method receives extra side information or tuning not available to baselines.
- Prefer BER/BLER curves across SNR, ablation studies, complexity/runtime analysis, and sensitivity to channel mismatch.

Practicality:
- Identify hidden overhead from pilots, side information, parameter search, feedback, synchronization, or channel estimation.
- Check whether complexity scales with subcarriers, antennas, users, paths, candidate parameters, or Monte Carlo trials.
- Distinguish simulation-only promise from implementable system design.

