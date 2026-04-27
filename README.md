# AFDM MATLAB Simulation Workspace

This workspace contains a MATLAB implementation of an AFDM transmit-channel-receive simulation chain. The earlier Python scaffold has been removed so the project is centered on the MATLAB workflow.

## Layout

- `matlab_afdm/afdm_config.m`: shared simulation, waveform, modulation, and channel configuration.
- `matlab_afdm/main_simulation.m`: single-run end-to-end simulation.
- `matlab_afdm/simulation_loop.m`: SNR sweep or Monte Carlo simulation entry point.
- `matlab_afdm/transmitter/`: random bit generation, QAM/PSK modulation, IDAFT modulation, CPP insertion, and PAPR calculation.
- `matlab_afdm/channel/`: multipath Doppler channel and AWGN.
- `matlab_afdm/receive/`: CPP removal, DAFT demodulation, effective-channel estimation, MMSE equalization, symbol decision, and BER counting.

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

Simulation mode and channel parameters are controlled in `matlab_afdm/afdm_config.m`.
