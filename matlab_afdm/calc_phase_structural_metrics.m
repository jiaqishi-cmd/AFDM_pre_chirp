function metrics = calc_phase_structural_metrics(phi_m, cfg)
%CALC_PHASE_STRUCTURAL_METRICS Compatibility wrapper.
%   Prefer afdm.metrics.phase_structural(phi_m, cfg) for new code.

    if nargin < 2
        cfg = struct();
    end
    metrics = afdm.metrics.phase_structural(phi_m, cfg);
end
