# AFDM MATLAB Simulation Workspace

This workspace contains a MATLAB implementation of an AFDM transmit-channel-receive simulation chain. The earlier Python scaffold has been removed so the project is centered on the MATLAB workflow.

## Layout

- `matlab_afdm/afdm_config.m`: shared simulation, waveform, modulation, and channel configuration.
- `matlab_afdm/main_simulation.m`: single-run end-to-end simulation.
- `matlab_afdm/simulation_loop.m`: SNR sweep or Monte Carlo simulation entry point.
- `matlab_afdm/+afdm/+tx/`: random bit generation, QAM/PSK modulation, IDAFT modulation, CPP insertion, and PAPR calculation.
- `matlab_afdm/+afdm/+channel/`: channel profile generation, multipath Doppler channel, and AWGN.
- `matlab_afdm/+afdm/+rx/`: CPP removal, DAFT demodulation, effective-channel estimation, equalization, symbol decision, and BER counting.
- `matlab_afdm/+afdm/+chirp/`: pre-chirp profiles, grouped c2 pattern builders, and frame-level profile selection.
- `matlab_afdm/+afdm/+analysis/`: reusable channel-matrix, distance, Phi, and Bemani-equation analysis helpers.
- `matlab_afdm/+afdm/+delta/`: structured, random, and recursive delta-set generators.
- `matlab_afdm/experiments/reproduction/`: formal reproduction and resume runners.
- Legacy `transmitter/`, `channel/`, `receive/`, and `pre_chirp/` wrapper directories were removed from the current branch after archiving older code on GitHub branches `legacy/pre-package-migration` and `legacy/pre-rx-channel-package`.
- Pre-chirp profile construction and frame-level selection now live in `matlab_afdm/+afdm/+chirp/`. GPS uses the paper-style candidate set, while the proposed profile keeps the original AFDM `c2` and greedily applies small group-wise perturbations from `{0, -delta, +delta}`.

## Run

From MATLAB:

```matlab
cd matlab_afdm
main_simulation
```

For SNR sweep or Monte Carlo evaluation:

```matlab
cd matlab_afdm
simulation_loop
```

For a baseline/GPS/proposed comparison:

```matlab
cd matlab_afdm
run_scheme_comparison
```

For proposed perturbation strength sweeps:

```matlab
cd matlab_afdm
run_delta_sweep
```

For higher-sample BER-only comparisons:

```matlab
cd matlab_afdm
run_ber_comparison
```

For Bemani two-path key-equation searches comparing uniform c2 and fixed
GPS c2,m group patterns:

```matlab
cd matlab_afdm
main_bemani_gps_key_equation_search
```

For stronger GPS-only Phi(delta) rank-loss searches and follow-up BER-SNR
validation:

```matlab
cd matlab_afdm
main_gps_unique_rank_loss_search
main_gps_bestcase_ber_snr
```

For adaptive BER runs with target-error stopping:

```matlab
cd matlab_afdm
run_adaptive_ber_comparison
```

Simulation mode and channel parameters are controlled in `matlab_afdm/afdm_config.m`.
By default, BER comparisons reuse the configured channel realization and pair
schemes with the same per-frame random seed. Set
`config.simulation.refresh_channel_per_frame = true` or pass
`options.refresh_channel_per_frame = true` to BER comparison helpers when you
want each simulated frame to draw a fresh channel profile.

Available channel profiles include the lightweight `random_3path` demo channel and `bemani_21path`, which follows the 21-path LTV setup used in Bemani et al. Fig. 5.

Example Bemani-style BER run:

```matlab
opts.channel_profile = 'bemani_21path';
opts.M_mod = 4;
opts.modType = 'qam';
run_adaptive_ber_comparison([15 20 25], [200 100 20], [5000 20000 50000], ...
    {'baseline'}, opts)
```
