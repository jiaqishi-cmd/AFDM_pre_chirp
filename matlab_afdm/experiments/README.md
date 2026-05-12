# MATLAB Experiment Scripts

This directory contains experiment entry scripts and plotting scripts. Core AFDM
functions remain in the `matlab_afdm` root and submodules so the simulation path
stays stable.

Folders:

- `bemani/`: Bemani key-equation, Phi(delta), and GPS rank-risk searches.
- `caseA/`: Case A fixed-channel stress tests and theta scans.
- `random_channel/`: Choi-style random doubly selective channel BER.
- `papr/`: PAPR CCDF, proposed delta sweep, and proposed BER/CCDF curves.
- `complexity/`: PAPR-search strategy and partial waveform reuse studies.
- `structural/`: Channel-independent c2 structural-risk metrics.

Run from `matlab_afdm` root, for example:

```matlab
setup_paths(pwd);
run('experiments/complexity/run_partial_reuse_theory_and_timing.m');
```

Each moved script uses `find_afdm_root` to locate the source root before calling
`setup_paths`, so running through `run('experiments/...')` keeps result output in
the repository-level `results/` directory.
