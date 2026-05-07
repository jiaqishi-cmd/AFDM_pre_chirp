function metrics = compute_bemani_equation_error(delta, d_gps, L, ldiff, N, c2_values)
%COMPUTE_BEMANI_EQUATION_ERROR Evaluate Bemani key equation errors.

    delta = delta(:);
    d_gps = d_gps(:);
    z = (0:N-1).';
    z_minus = mod(z - L, N) + 1;
    z_plus = mod(z + L, N) + 1;

    ratio_delta = delta.^2 ./ (delta(z_minus) .* delta(z_plus));
    path_phase = exp(1i * 2 * pi / N * ldiff * L);

    metrics.ratio_delta = ratio_delta;
    metrics.path_phase = path_phase;

    names = fieldnames(c2_values);
    for idx = 1:numel(names)
        name = names{idx};
        phase_uniform = path_phase * exp(1i * 4 * pi * c2_values.(name) * L^2);
        errors = abs(ratio_delta - phase_uniform);
        metrics.(['min_E_' name]) = min(errors);
        metrics.(['mean_E_' name]) = mean(errors);
        metrics.(['max_E_' name]) = max(errors);
    end

    chi_gps = d_gps.^2 ./ (d_gps(z_minus) .* d_gps(z_plus));
    phase_gps = path_phase * chi_gps;
    errors_gps = abs(ratio_delta - phase_gps);

    metrics.chi_gps = chi_gps;
    metrics.phase_gps = phase_gps;
    metrics.min_E_GPS = min(errors_gps);
    metrics.mean_E_GPS = mean(errors_gps);
    metrics.max_E_GPS = max(errors_gps);
end
