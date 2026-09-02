%CHECK_ENVIRONMENT Verify that the AFDM simulation workspace can run locally.

rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(rootDir);
setup_paths(rootDir);

fprintf('========== AFDM workspace environment check ==========\n');
fprintf('Workspace root: %s\n', rootDir);
fprintf('MATLAB version: %s\n', version);

requiredFunctions = { ...
    'afdm_config', ...
    'simulate_frame', ...
    'afdm.tx.engine', ...
    'afdm.channel.multipath', ...
    'afdm.channel.add_awgn', ...
    'afdm.rx.engine', ...
    'qammod'};

for idx = 1:numel(requiredFunctions)
    name = requiredFunctions{idx};
    location = which(name);
    if isempty(location)
        error('Required function not found: %s', name);
    end
    fprintf('OK: %-28s -> %s\n', name, location);
end

cfg = afdm_config();
cfg.waveform.NumSubcarriers = 16;
cfg.modulation.M_mod = 4;
cfg.simulation.num_frames = 1;

frame = simulate_frame(cfg, cfg.simulation.random_seed);
fprintf('Smoke test: PAPR = %.2f dB, errors = %d / %d, BER = %.3g\n', ...
    frame.papr, frame.err_bits, frame.total_bits, frame.ber);

[signalCpp, ~, txBits] = afdm.tx.engine(cfg);
rSignal = afdm.channel.multipath(signalCpp, cfg);
[rSignal, noiseVar] = afdm.channel.add_awgn(rSignal, cfg);
cfg.channel.noise_var = noiseVar;
[~, errBits, totalBits] = afdm_rx_engine(rSignal, cfg, txBits);
fprintf('Legacy wrapper: errors = %d / %d\n', errBits, totalBits);

fprintf('Environment check complete.\n');
