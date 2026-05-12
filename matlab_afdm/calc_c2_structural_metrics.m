function metrics = calc_c2_structural_metrics(c2_vec, cfg)
%CALC_C2_STRUCTURAL_METRICS Compatibility wrapper.
%   Prefer afdm.metrics.c2_structural(c2_vec, cfg) for new code.

    if nargin < 2
        cfg = struct();
    end
    metrics = afdm.metrics.c2_structural(c2_vec, cfg);
end
