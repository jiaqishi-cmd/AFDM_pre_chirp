function corrMetrics = c2_candidate_correlation(c2VecList, cfg)
%C2_CANDIDATE_CORRELATION Compute similarity between c2 candidate patterns.
%   rho_ab = abs(mean(exp(1j*2*pi*(c2_a-c2_b).*m.^2))).

    if nargin < 2
        cfg = struct();
    end

    C = normalize_candidate_input(c2VecList);
    M = size(C, 1);
    K = size(C, 2);
    if isstruct(cfg) && isfield(cfg, 'M') && ~isempty(cfg.M) && cfg.M ~= M
        error('cfg.M does not match c2VecList length.');
    end
    if isstruct(cfg) && isfield(cfg, 'N') && ~isempty(cfg.N) && cfg.N ~= M
        error('cfg.N does not match c2VecList length.');
    end

    m = (0:M-1).';
    rho = eye(K);
    for a = 1:K
        for b = a+1:K
            phaseDiff = exp(1i * 2 * pi * (C(:, a) - C(:, b)) .* (m.^2));
            rho(a, b) = abs(mean(phaseDiff));
            rho(b, a) = rho(a, b);
        end
    end

    offdiag = ~eye(K);
    if K > 1
        maxRho = max(rho(offdiag));
        meanRho = mean(rho(offdiag));
    else
        maxRho = NaN;
        meanRho = NaN;
    end

    corrMetrics = struct();
    corrMetrics.rho_matrix = rho;
    corrMetrics.max_rho_offdiag = maxRho;
    corrMetrics.mean_rho_offdiag = meanRho;
    corrMetrics.min_separation = 1 - maxRho;
end

function C = normalize_candidate_input(c2VecList)
    if iscell(c2VecList)
        K = numel(c2VecList);
        C = zeros(numel(c2VecList{1}), K);
        for idx = 1:K
            C(:, idx) = c2VecList{idx}(:);
        end
    else
        C = c2VecList;
        if isvector(C)
            C = C(:);
        end
    end
end
