function metrics = c2_structural(c2Vec, cfg)
%C2_STRUCTURAL Compute channel-independent structural-risk metrics for c2.
%   metrics = afdm.metrics.c2_structural(c2Vec, cfg)
%
%   R_struct is a channel-independent structural-risk indicator. It does not
%   prove full-diversity loss directly; it only flags phase patterns that may
%   align more easily with constellation difference structures.

    if nargin < 2
        cfg = struct();
    end

    c2Vec = c2Vec(:);
    M = infer_length(c2Vec, cfg);
    if numel(c2Vec) ~= M
        error('c2Vec length must match cfg.M or cfg.N.');
    end

    m = (0:M-1).';
    phi = exp(1i * 2 * pi .* c2Vec .* (m.^2));
    metrics = afdm.metrics.phase_structural(phi, cfg);
    metrics.c2_vec = c2Vec;
end

function M = infer_length(c2Vec, cfg)
    if isstruct(cfg) && isfield(cfg, 'M') && ~isempty(cfg.M)
        M = cfg.M;
    elseif isstruct(cfg) && isfield(cfg, 'N') && ~isempty(cfg.N)
        M = cfg.N;
    else
        M = numel(c2Vec);
    end
end
