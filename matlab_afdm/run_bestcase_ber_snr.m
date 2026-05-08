function results = run_bestcase_ber_snr(best_case, options)
%RUN_BESTCASE_BER_SNR 对搜索得到的 best_case 做 BER-SNR 验证。
%   使用固定两径支撑，路径增益每帧 Rayleigh 随机。当前工程没有 ML
%   检测器，因此默认使用现有 MMSE 接收链路。

    if nargin < 2
        options = struct();
    end

    snrValues = get_option(options, 'snr_values', 0:2:30);
    numFrames = get_option(options, 'num_frames', 500);
    rngSeed = get_option(options, 'seed', 20260427);

    schemes = {'baseline', 'GPS', 'proposed'};
    N = best_case.N;
    V = best_case.V;
    c1 = best_case.c1;
    c2_base = best_case.c2_base;
    gpsPattern = parse_pattern(best_case.gps_pattern);
    propPattern = parse_pattern(best_case.proposed_pattern);
    proposed_delta = c2_base / 16;
    [c2_gps, ~] = build_c2m_gps_pattern(N, V, gpsPattern);
    [c2_prop, ~] = build_c2m_proposed_pattern(N, V, propPattern, c2_base, proposed_delta);

    results.best_case = best_case;
    results.snr_values = snrValues;
    results.num_frames = numFrames;
    results.schemes = schemes;
    results.ber = zeros(numel(snrValues), numel(schemes));
    results.error_bits = zeros(numel(snrValues), numel(schemes));
    results.total_bits = zeros(numel(snrValues), numel(schemes));

    fprintf('========== BER-SNR on GPS best case ==========\n');
    fprintf('N=%d V=%d L=%d l2=%d alpha2=%d delta=%s\n', ...
        N, V, best_case.L, best_case.l2, best_case.alpha2, best_case.delta_type);

    for snrIdx = 1:numel(snrValues)
        for schemeIdx = 1:numel(schemes)
            totalErr = 0;
            totalBits = 0;
            for frameIdx = 1:numFrames
                cfg = build_frame_config(best_case, schemes{schemeIdx}, c2_base, c2_gps, c2_prop, c1, snrValues(snrIdx));
                rng(rngSeed + 100000 * snrIdx + 1000 * schemeIdx + frameIdx, 'twister');
                cfg.channel.chan_coef = (randn(1, 2) + 1i * randn(1, 2)) / sqrt(4);
                frame = simulate_frame(cfg, rngSeed + 1000000 * snrIdx + 10000 * schemeIdx + frameIdx);
                totalErr = totalErr + frame.err_bits;
                totalBits = totalBits + frame.total_bits;
            end
            results.error_bits(snrIdx, schemeIdx) = totalErr;
            results.total_bits(snrIdx, schemeIdx) = totalBits;
            results.ber(snrIdx, schemeIdx) = totalErr / totalBits;
            fprintf('%-8s | SNR %5.1f dB | errors %6d/%-8d | BER %.3e\n', ...
                schemes{schemeIdx}, snrValues(snrIdx), totalErr, totalBits, results.ber(snrIdx, schemeIdx));
        end
    end

    results.diversity_order = estimate_diversity_order(results.snr_values, results.ber);
    plot_ber_results(results);
end

function cfg = build_frame_config(best_case, scheme, c2_base, c2_gps, c2_prop, c1, snrDb)
    cfg = afdm_config();
    cfg.waveform.NumSubcarriers = best_case.N;
    cfg.waveform.CPPLength = max(best_case.l2, 1);
    cfg.waveform.c1 = c1;
    cfg.modulation.M_mod = 2;
    cfg.modulation.modType = 'psk';
    cfg.channel.multipath = true;
    cfg.channel.add_noise = true;
    cfg.channel.snr_db = snrDb;
    cfg.channel.delay_taps = [best_case.l1, best_case.l2];
    cfg.channel.doppler_taps = [best_case.alpha1, best_case.alpha2];
    cfg.channel.doppler_freq = cfg.channel.doppler_taps / cfg.waveform.NumSubcarriers;

    switch lower(scheme)
        case 'baseline'
            c2 = c2_base;
        case 'gps'
            c2 = c2_gps;
        case 'proposed'
            c2 = c2_prop;
        otherwise
            error('Unknown scheme: %s', scheme);
    end

    cfg.pre_chirp.scheme = 'baseline';
    cfg.pre_chirp.profile.scheme = 'baseline';
    cfg.pre_chirp.profile.c2 = c2;
    cfg.waveform.c2 = c2;
end

function div = estimate_diversity_order(snrDb, ber)
    snrLinear = 10 .^ (snrDb(:) / 10);
    div = zeros(1, size(ber, 2));
    for schemeIdx = 1:size(ber, 2)
        idx = snrDb(:) >= 15 & ber(:, schemeIdx) > 0;
        if nnz(idx) >= 2
            p = polyfit(log10(snrLinear(idx)), log10(ber(idx, schemeIdx)), 1);
            div(schemeIdx) = -p(1);
        else
            div(schemeIdx) = NaN;
        end
    end
end

function plot_ber_results(results)
    outputDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    figure('Name', 'GPS best-case BER-SNR', 'Color', 'w');
    semilogy(results.snr_values, results.ber, 'o-', 'LineWidth', 2);
    grid on;
    xlabel('SNR (dB)');
    ylabel('BER');
    title('BER-SNR on GPS near-rank-loss best case');
    legend(results.schemes, 'Location', 'southwest');
    saveas(gcf, fullfile(outputDir, ['gps_bestcase_ber_snr_' timestamp '.png']));
    save(fullfile(outputDir, 'gps_bestcase_ber_snr.mat'), 'results');
    fprintf('Diversity estimates: baseline %.3g | GPS %.3g | proposed %.3g\n', ...
        results.diversity_order(1), results.diversity_order(2), results.diversity_order(3));
end

function pattern = parse_pattern(text)
    if isnumeric(text)
        pattern = text;
        return;
    end
    nums = regexp(char(text), '\d+', 'match');
    pattern = str2double(nums);
end

function value = get_option(options, name, defaultValue)
    if isstruct(options) && isfield(options, name)
        value = options.(name);
    else
        value = defaultValue;
    end
end
