function metrics = calc_phase_structural_metrics(phi_m, cfg)
%CALC_PHASE_STRUCTURAL_METRICS 直接基于 phase mask 计算结构风险指标。
%   该函数不使用信道 h、delay、Doppler，只看相位分布和星座差分相位对齐。

    if nargin < 2
        cfg = struct();
    end

    phi_m = phi_m(:);
    M = numel(phi_m);
    numBins = get_option(cfg, 'num_bins', 64);
    phaseTol = get_option(cfg, 'phase_tol', pi / 24);
    modType = upper(get_option(cfg, 'modType', 'BPSK'));
    weights = get_option(cfg, 'struct_weights', [0.3, 0.2, 0.5]);
    weights = weights(:).' / sum(weights);

    phase = mod(angle(phi_m), 2 * pi);
    edges = linspace(0, 2 * pi, numBins + 1);
    counts = histcounts(phase, edges, 'Normalization', 'probability');
    p = counts(:);

    phaseEntropy = -sum(p .* log(p + eps)) / log(numBins);
    phaseRisk = 1 - phaseEntropy;

    effBins = sum(counts > 0);
    effPhaseRatio = effBins / numBins;
    uniquePhaseRisk = 1 - effPhaseRatio;

    targetAngles = get_target_angles(modType);
    phaseDiff = angle(phi_m * phi_m');
    nontrivial = ~eye(M);
    minDist = inf(M, M);
    for idx = 1:numel(targetAngles)
        minDist = min(minDist, abs(local_wrap_to_pi(phaseDiff - targetAngles(idx))));
    end
    aligned = minDist < phaseTol;
    alignmentRatio = sum(aligned(nontrivial), 'all') / nnz(nontrivial);

    Rstruct = weights(1) * phaseRisk + ...
        weights(2) * uniquePhaseRisk + ...
        weights(3) * alignmentRatio;

    metrics = struct();
    metrics.phase_entropy = phaseEntropy;
    metrics.phase_degeneracy_risk = phaseRisk;
    metrics.eff_phase_bins = effBins;
    metrics.eff_phase_ratio = effPhaseRatio;
    metrics.unique_phase_risk = uniquePhaseRisk;
    metrics.constellation_alignment_ratio = alignmentRatio;
    metrics.constellation_alignment_risk = alignmentRatio;
    metrics.R_struct = Rstruct;
    metrics.phase = phase;
    metrics.phase_hist_counts = counts;
    metrics.phase_hist_edges = edges;
    metrics.num_bins = numBins;
    metrics.phase_tol = phaseTol;
    metrics.modType = modType;
    metrics.struct_weights = weights;
end

function targetAngles = get_target_angles(modType)
    switch upper(modType)
        case 'BPSK'
            targetAngles = [0, pi];
        case {'QPSK', '16QAM'}
            % TODO: 16QAM 可扩展为真实 constellation difference angle set。
            targetAngles = [0, pi/2, pi, 3*pi/2];
        otherwise
            error('Unsupported modType: %s', modType);
    end
end

function y = local_wrap_to_pi(x)
    y = mod(x + pi, 2 * pi) - pi;
end

function value = get_option(cfg, name, defaultValue)
    if isstruct(cfg) && isfield(cfg, name) && ~isempty(cfg.(name))
        value = cfg.(name);
    else
        value = defaultValue;
    end
end
