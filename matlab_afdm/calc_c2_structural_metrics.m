function metrics = calc_c2_structural_metrics(c2_vec, cfg)
%CALC_C2_STRUCTURAL_METRICS 计算 c2 pattern 的信道无关结构风险指标。
%   metrics = calc_c2_structural_metrics(c2_vec, cfg)
%
%   R_struct is a channel-independent structural-risk indicator.
%   It does not prove full-diversity loss directly.
%   It is used to identify phase patterns that may be more likely to align
%   with constellation difference structures.

    if nargin < 2
        cfg = struct();
    end

    c2_vec = c2_vec(:);
    M = infer_length(c2_vec, cfg);
    if numel(c2_vec) ~= M
        error('c2_vec length must match cfg.M or cfg.N.');
    end

    m = (0:M-1).';
    phi_m = exp(1i * 2 * pi .* c2_vec .* (m.^2));
    metrics = calc_phase_structural_metrics(phi_m, cfg);
    metrics.c2_vec = c2_vec;
end

function M = infer_length(c2_vec, cfg)
    if isstruct(cfg) && isfield(cfg, 'M') && ~isempty(cfg.M)
        M = cfg.M;
    elseif isstruct(cfg) && isfield(cfg, 'N') && ~isempty(cfg.N)
        M = cfg.N;
    else
        M = numel(c2_vec);
    end
end
