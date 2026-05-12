function corr_metrics = calc_c2_candidate_correlation(c2_vec_list, cfg)
%CALC_C2_CANDIDATE_CORRELATION Compatibility wrapper.
%   Prefer afdm.metrics.c2_candidate_correlation(c2_vec_list, cfg) for new code.

    if nargin < 2
        cfg = struct();
    end
    corr_metrics = afdm.metrics.c2_candidate_correlation(c2_vec_list, cfg);
end
