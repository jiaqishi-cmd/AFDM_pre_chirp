function corr_metrics = calc_c2_candidate_correlation(c2_vec_list, cfg)
%CALC_C2_CANDIDATE_CORRELATION 计算多个 c2 candidate pattern 之间的相似度。
%   rho_ab = abs(mean(exp(1j*2*pi*(c2_a-c2_b).*m.^2))).

    if nargin < 2
        cfg = struct();
    end

    C = normalize_candidate_input(c2_vec_list);
    M = size(C, 1);
    K = size(C, 2);
    if isstruct(cfg) && isfield(cfg, 'M') && ~isempty(cfg.M) && cfg.M ~= M
        error('cfg.M does not match c2_vec_list length.');
    end
    if isstruct(cfg) && isfield(cfg, 'N') && ~isempty(cfg.N) && cfg.N ~= M
        error('cfg.N does not match c2_vec_list length.');
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

    corr_metrics = struct();
    corr_metrics.rho_matrix = rho;
    corr_metrics.max_rho_offdiag = maxRho;
    corr_metrics.mean_rho_offdiag = meanRho;
    corr_metrics.min_separation = 1 - maxRho;
end

function C = normalize_candidate_input(c2_vec_list)
    if iscell(c2_vec_list)
        K = numel(c2_vec_list);
        C = zeros(numel(c2_vec_list{1}), K);
        for idx = 1:K
            C(:, idx) = c2_vec_list{idx}(:);
        end
    else
        C = c2_vec_list;
        if isvector(C)
            C = C(:);
        end
    end
end
