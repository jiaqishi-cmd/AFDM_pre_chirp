function metrics = phase_structural(phi, cfg)
%PHASE_STRUCTURAL Compute structural-risk metrics from a phase mask.
%   This function uses only the phase distribution and constellation
%   difference directions. It does not use channel gains, delay, or Doppler.

    if nargin < 2
        cfg = struct();
    end

    phi = phi(:);
    M = numel(phi);
    numBins = get_option(cfg, 'num_bins', 64);
    phaseTol = get_option(cfg, 'phase_tol', pi / 24);
    modType = upper(get_option(cfg, 'modType', 'BPSK'));
    weights = get_option(cfg, 'struct_weights', [0.3, 0.2, 0.5]);
    weights = weights(:).' / sum(weights);

    phase = mod(angle(phi), 2 * pi);
    edges = linspace(0, 2 * pi, numBins + 1);
    counts = histcounts(phase, edges, 'Normalization', 'probability');
    p = counts(:);

    phaseEntropy = -sum(p .* log(p + eps)) / log(numBins);
    phaseRisk = 1 - phaseEntropy;

    effBins = sum(counts > 0);
    effPhaseRatio = effBins / numBins;
    uniquePhaseRisk = 1 - effPhaseRatio;

    targetAngles = target_angles(modType);
    phaseDiff = angle(phi * phi');
    nontrivial = ~eye(M);
    minDist = inf(M, M);
    for idx = 1:numel(targetAngles)
        minDist = min(minDist, abs(wrap_to_pi(phaseDiff - targetAngles(idx))));
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

function targetAngles = target_angles(modType)
    switch upper(modType)
        case 'BPSK'
            targetAngles = [0, pi];
        case {'QPSK', '16QAM'}
            % TODO: extend 16QAM to its exact constellation-difference set.
            targetAngles = [0, pi/2, pi, 3*pi/2];
        otherwise
            error('Unsupported modType: %s', modType);
    end
end

function y = wrap_to_pi(x)
    y = mod(x + pi, 2 * pi) - pi;
end

function value = get_option(cfg, name, defaultValue)
    if isstruct(cfg) && isfield(cfg, name) && ~isempty(cfg.(name))
        value = cfg.(name);
    else
        value = defaultValue;
    end
end
