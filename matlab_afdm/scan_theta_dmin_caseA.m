function scan = scan_theta_dmin_caseA(theta_list, delta_set, options)
%SCAN_THETA_DMIN_CASEA 对 Case A 扫描两径增益相位并计算有限星座 dmin。
%   主指标是 dmin2，不用 total H_eff rank 作为结论；rank/cond 仅辅助排查。

    if nargin < 3
        options = struct();
    end

    N = get_option(options, 'N', 64);
    V = get_option(options, 'V', 4);
    M = N / V; %#ok<NASGU>
    alphaMax = get_option(options, 'alpha_max', 3);
    c1 = get_option(options, 'c1', (2 * alphaMax + 1) / (2 * N));
    c2Base = get_option(options, 'c2_base', sqrt(2) / (10 * N));
    gpsPattern = get_option(options, 'gps_pattern', [1 2 1 2]);
    proposedPattern = get_option(options, 'proposed_pattern', gpsPattern);
    proposedDelta = get_option(options, 'proposed_delta', c2Base / 16);

    l1 = 0;
    alpha1 = 0;
    l2 = 2;
    alpha2 = 2;

    [c2Gps, ~] = build_c2m_gps_pattern(N, V, gpsPattern);
    hasProposed = exist('build_c2m_proposed_pattern', 'file') == 2;
    if hasProposed
        [c2Prop, ~] = build_c2m_proposed_pattern(N, V, proposedPattern, c2Base, proposedDelta);
    else
        c2Prop = NaN;
    end

    paths.base.H1 = build_H_path_general_c2m(N, c1, c2Base, l1, alpha1);
    paths.base.H2 = build_H_path_general_c2m(N, c1, c2Base, l2, alpha2);
    paths.gps.H1 = build_H_path_general_c2m(N, c1, c2Gps, l1, alpha1);
    paths.gps.H2 = build_H_path_general_c2m(N, c1, c2Gps, l2, alpha2);
    if hasProposed
        paths.prop.H1 = build_H_path_general_c2m(N, c1, c2Prop, l1, alpha1);
        paths.prop.H2 = build_H_path_general_c2m(N, c1, c2Prop, l2, alpha2);
    end

    pre.base.Y1 = paths.base.H1 * delta_set.delta;
    pre.base.Y2 = paths.base.H2 * delta_set.delta;
    pre.gps.Y1 = paths.gps.H1 * delta_set.delta;
    pre.gps.Y2 = paths.gps.H2 * delta_set.delta;
    if hasProposed
        pre.prop.Y1 = paths.prop.H1 * delta_set.delta;
        pre.prop.Y2 = paths.prop.H2 * delta_set.delta;
    end

    numTheta = numel(theta_list);
    rows = repmat(make_empty_row(), numTheta, 1);
    bestDeltaGps = [];

    for thetaIdx = 1:numTheta
        theta = theta_list(thetaIdx);
        h1 = 1 / sqrt(2);
        h2 = exp(1i * theta) / sqrt(2);

        base = dmin_from_precomputed(pre.base.Y1, pre.base.Y2, h1, h2, delta_set);
        gps = dmin_from_precomputed(pre.gps.Y1, pre.gps.Y2, h1, h2, delta_set);
        if hasProposed
            prop = dmin_from_precomputed(pre.prop.Y1, pre.prop.Y2, h1, h2, delta_set);
        else
            prop = make_nan_dmin();
        end

        Hbase = h1 * paths.base.H1 + h2 * paths.base.H2;
        Hgps = h1 * paths.gps.H1 + h2 * paths.gps.H2;
        if hasProposed
            Hprop = h1 * paths.prop.H1 + h2 * paths.prop.H2;
        else
            Hprop = NaN(N);
        end

        row = make_empty_row();
        row.theta = theta;
        row.theta_over_pi = theta / pi;
        row.dmin2_base = base.dmin2;
        row.dmin2_GPS = gps.dmin2;
        row.dmin2_prop = prop.dmin2;
        row.delta_type_base = base.best_delta_type;
        row.delta_type_GPS = gps.best_delta_type;
        row.delta_type_prop = prop.best_delta_type;
        row.delta_index_base = base.best_delta_index;
        row.delta_index_GPS = gps.best_delta_index;
        row.delta_index_prop = prop.best_delta_index;
        row.ratio_GPS_base = gps.dmin2 / (base.dmin2 + eps);
        row.ratio_GPS_prop = gps.dmin2 / (prop.dmin2 + eps);
        row.gap_base_GPS_dB = 10 * log10((base.dmin2 + eps) / (gps.dmin2 + eps));
        row.gap_prop_GPS_dB = 10 * log10((prop.dmin2 + eps) / (gps.dmin2 + eps));
        row.rank_H_base = rank(Hbase, 1e-10);
        row.rank_H_GPS = rank(Hgps, 1e-10);
        row.rank_H_prop = rank(Hprop, 1e-10);
        row.cond_H_base = cond_from_svd(Hbase);
        row.cond_H_GPS = cond_from_svd(Hgps);
        row.cond_H_prop = cond_from_svd(Hprop);
        row.score = (base.dmin2 + prop.dmin2) / (2 * (gps.dmin2 + eps));
        rows(thetaIdx) = row;

        if isempty(bestDeltaGps) || row.score > rows(bestDeltaGps.row_index).score
            bestDeltaGps.row_index = thetaIdx;
            bestDeltaGps.delta = gps.best_delta;
        end
    end

    phase_scan_table = struct2table(rows);
    best = select_best_theta(phase_scan_table);
    best_delta_GPS = delta_set.delta(:, phase_scan_table.delta_index_GPS(best.row_index)); %#ok<NASGU>

    scan.phase_scan_table = phase_scan_table;
    scan.best_case = best;
    scan.best_delta_GPS = best_delta_GPS;
    scan.theta_list = theta_list;
    scan.gps_pattern = gpsPattern;
    scan.proposed_pattern = proposedPattern;
    scan.has_proposed = hasProposed;
    scan.paths = paths;
end

function out = dmin_from_precomputed(Y1, Y2, h1, h2, delta_set)
    Y = h1 * Y1 + h2 * Y2;
    metric = sum(abs(Y).^2, 1) ./ delta_set.norm2;
    [out.dmin2, bestIdx] = min(real(metric));
    out.best_delta = delta_set.delta(:, bestIdx);
    out.best_delta_type = delta_set.type(bestIdx);
    out.best_delta_index = bestIdx;
end

function best = select_best_theta(T)
    thresholdBase = 0.1;
    thresholdProp = 0.1;
    valid = T.dmin2_base > thresholdBase & T.dmin2_prop > thresholdProp;
    if ~any(valid)
        thresholdBase = 0.03;
        thresholdProp = 0.03;
        valid = T.dmin2_base > thresholdBase & T.dmin2_prop > thresholdProp;
    end
    if ~any(valid)
        valid = true(height(T), 1);
    end

    score = T.score;
    score(~valid) = -Inf;
    [~, idx] = max(score);

    best = struct();
    best.row_index = idx;
    best.theta = T.theta(idx);
    best.theta_over_pi = T.theta_over_pi(idx);
    best.dmin2_base = T.dmin2_base(idx);
    best.dmin2_GPS = T.dmin2_GPS(idx);
    best.dmin2_prop = T.dmin2_prop(idx);
    best.delta_type_GPS = char(T.delta_type_GPS(idx));
    best.delta_index_GPS = T.delta_index_GPS(idx);
    best.gap_base_GPS_dB = T.gap_base_GPS_dB(idx);
    best.gap_prop_GPS_dB = T.gap_prop_GPS_dB(idx);
    best.score = T.score(idx);
    best.threshold_base = thresholdBase;
    best.threshold_prop = thresholdProp;
end

function row = make_empty_row()
    row = struct( ...
        'theta', 0, ...
        'theta_over_pi', 0, ...
        'dmin2_base', NaN, ...
        'dmin2_GPS', NaN, ...
        'dmin2_prop', NaN, ...
        'delta_type_base', "", ...
        'delta_type_GPS', "", ...
        'delta_type_prop', "", ...
        'delta_index_base', NaN, ...
        'delta_index_GPS', NaN, ...
        'delta_index_prop', NaN, ...
        'ratio_GPS_base', NaN, ...
        'ratio_GPS_prop', NaN, ...
        'gap_base_GPS_dB', NaN, ...
        'gap_prop_GPS_dB', NaN, ...
        'rank_H_base', NaN, ...
        'rank_H_GPS', NaN, ...
        'rank_H_prop', NaN, ...
        'cond_H_base', NaN, ...
        'cond_H_GPS', NaN, ...
        'cond_H_prop', NaN, ...
        'score', NaN);
end

function dmin = make_nan_dmin()
    dmin.dmin2 = NaN;
    dmin.best_delta = NaN;
    dmin.best_delta_type = "";
    dmin.best_delta_index = NaN;
end

function c = cond_from_svd(H)
    s = svd(H);
    c = max(s) / max(min(s), eps);
end

function value = get_option(options, name, defaultValue)
    if isstruct(options) && isfield(options, name)
        value = options.(name);
    else
        value = defaultValue;
    end
end
