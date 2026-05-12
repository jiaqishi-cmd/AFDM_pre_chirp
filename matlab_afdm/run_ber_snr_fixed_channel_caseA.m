function ber_results = run_ber_snr_fixed_channel_caseA(theta, gps_pattern, options)
%RUN_BER_SNR_FIXED_CHANNEL_CASEA 固定 Case A 和两径相位跑 BER-SNR。

    if nargin < 3
        options = struct();
    end

    N = 64;
    V = 4;
    alphaMax = 3;
    c1 = (2 * alphaMax + 1) / (2 * N);

    best_case = struct();
    best_case.N = N;
    best_case.V = V;
    best_case.M = N / V;
    best_case.alpha_max = alphaMax;
    best_case.c1 = c1;
    best_case.c2_base = sqrt(2) / (10 * N);
    best_case.l1 = 0;
    best_case.alpha1 = 0;
    best_case.l2 = 2;
    best_case.alpha2 = 2;
    best_case.L = 16;
    best_case.delta_type = 'phase_scan';
    best_case.gps_pattern = gps_pattern;
    best_case.proposed_pattern = gps_pattern;

    snrValues = get_option(options, 'snr_values', 0:2:34);
    frames = get_option(options, 'num_frames', 2000);
    modulationList = get_option(options, 'modulation_list', [2 4]);
    seed = get_option(options, 'seed', 20260508);
    labelPrefix = get_option(options, 'label_prefix', 'caseA_besttheta');

    ber_results = struct();
    for idx = 1:numel(modulationList)
        MMod = modulationList(idx);
        opts = struct();
        opts.snr_values = snrValues;
        opts.num_frames = frames;
        opts.seed = seed + 1000 * idx;
        opts.gain_mode = 'fixed';
        opts.fixed_theta = theta;
        opts.M_mod = MMod;
        opts.label = sprintf('%s_M%d_theta_%0.4fpi', labelPrefix, MMod, theta / pi);

        result = run_bestcase_ber_snr(best_case, opts);
        result.unstable_points = result.error_bits < 20 & result.error_bits > 0;
        result.slope_reliable = estimate_reliable_slope(result);

        if MMod == 2
            ber_results.bpsk = result;
        elseif MMod == 4
            ber_results.qpsk = result;
        else
            ber_results.(sprintf('M%d', MMod)) = result;
        end
    end
end

function slope = estimate_reliable_slope(result)
    snrDb = result.snr_values(:);
    snrLinear = 10 .^ (snrDb / 10);
    slope = NaN(1, size(result.ber, 2));
    for schemeIdx = 1:size(result.ber, 2)
        err = result.error_bits(:, schemeIdx);
        ber = result.ber(:, schemeIdx);
        valid = snrDb >= 16 & ber > 0 & err >= 20;
        if nnz(valid) < 3
            valid = snrDb >= 12 & ber > 0;
        end
        if nnz(valid) >= 3
            p = polyfit(log10(snrLinear(valid)), log10(ber(valid)), 1);
            slope(schemeIdx) = -p(1);
        end
    end
end

function value = get_option(options, name, defaultValue)
    if isstruct(options) && isfield(options, name)
        value = options.(name);
    else
        value = defaultValue;
    end
end
