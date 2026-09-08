% TEST_PAIRWISE_DIVERSITY_METRICS
rootDir = find_afdm_root(fileparts(mfilename('fullpath')));
addpath(rootDir);
setup_paths(rootDir);

pathMatrices = {eye(4), circshift(eye(4), 1, 1)};
deltaX = [1; 0; 0; 0];
psi = afdm.analysis.build_pairwise_path_response(pathMatrices, deltaX);
metrics = afdm.analysis.pairwise_diversity_metrics(psi);

assert(isequal(size(psi), [4, 2]));
assert(metrics.rank == 2);
assert(metrics.min_singular_value > 0);

dependentPaths = {eye(4), 2 * eye(4)};
psiDependent = afdm.analysis.build_pairwise_path_response( ...
    dependentPaths, deltaX);
dependentMetrics = afdm.analysis.pairwise_diversity_metrics(psiDependent);
assert(dependentMetrics.rank == 1);
assert(dependentMetrics.gram_determinant < 1e-12);

deltaSet = afdm.delta.generate_sparse_difference_set(4, 2, 'psk', 4);
assert(size(deltaSet.delta, 2) == 3^4 - 1);
assert(all(deltaSet.norm2 > 0));

N = 8;
c1 = 5 / (2 * N);
c2 = afdm.chirp.build_gps_pattern(N, 4, [2 2 1 1]);
gains = [0.7 + 0.2i, -0.4 + 0.3i];
delays = [0, 2];
dopplers = [0, 2];
symbols = pskmod(mod((0:N-1).', 2), 2, pi / 2);

cfg = afdm_config();
cfg.waveform.NumSubcarriers = N;
cfg.waveform.CPPLength = max(delays);
cfg.waveform.c1 = c1;
cfg.waveform.c2 = c2;
cfg.channel.delay_taps = delays;
cfg.channel.doppler_taps = dopplers;
cfg.channel.doppler_freq = dopplers / N;
cfg.channel.chan_coef = gains;

tx = afdm.tx.idaft_mod(symbols, N, c1, c2);
txCpp = afdm.tx.add_cpp(tx, cfg.waveform.CPPLength, c1);
rxCpp = afdm.channel.multipath(txCpp, cfg);
rx = afdm.rx.remove_cpp(rxCpp, cfg.waveform.CPPLength);
yDaft = afdm.rx.daft_demod(rx, N, c1, c2);
hEff = afdm.rx.estimate_effective_channel( ...
    N, c1, c2, gains, delays, dopplers / N);
relativeError = norm(yDaft - hEff * symbols) / norm(yDaft);
assert(relativeError < 1e-10);

cycleOperator = circshift(eye(4), 1, 1);
cycleAnalysis = afdm.analysis.monomial_cycle_analysis( ...
    cycleOperator, [-2i; 2i]);
assert(cycleAnalysis.constant_shift == 1);
assert(cycleAnalysis.min_compatible_weight == 4);
assert(cycleAnalysis.num_compatible_eigenvectors >= 1);

[gpsRational, ~] = afdm.chirp.build_gps_pattern_variant( ...
    8, 4, [2 2 1 1], 'rational');
[gpsIrrational, ~] = afdm.chirp.build_gps_pattern_variant( ...
    8, 4, [2 2 1 1], 'irrational');
[gpsLegacy, ~] = afdm.chirp.build_gps_pattern(8, 4, [2 2 1 1]);
assert(norm(gpsRational - gpsLegacy) < 1e-12);
assert(norm(gpsIrrational - gpsRational) > 0);

fprintf('test_pairwise_diversity_metrics passed.\n');
