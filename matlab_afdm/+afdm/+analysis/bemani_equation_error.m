function metrics = bemani_equation_error(delta, dGps, L, delayDiff, N, c2Values)
%BEMANI_EQUATION_ERROR Evaluate Bemani key-equation errors.

    delta = delta(:);
    dGps = dGps(:);
    z = (0:N-1).';
    zMinus = mod(z - L, N) + 1;
    zPlus = mod(z + L, N) + 1;

    ratioDelta = delta.^2 ./ (delta(zMinus) .* delta(zPlus));
    pathPhase = exp(1i * 2 * pi / N * delayDiff * L);

    metrics.ratio_delta = ratioDelta;
    metrics.path_phase = pathPhase;

    names = fieldnames(c2Values);
    for idx = 1:numel(names)
        name = names{idx};
        phaseUniform = pathPhase * exp(1i * 4 * pi * c2Values.(name) * L^2);
        errors = abs(ratioDelta - phaseUniform);
        metrics.(['min_E_' name]) = min(errors);
        metrics.(['mean_E_' name]) = mean(errors);
        metrics.(['max_E_' name]) = max(errors);
    end

    chiGps = dGps.^2 ./ (dGps(zMinus) .* dGps(zPlus));
    phaseGps = pathPhase * chiGps;
    errorsGps = abs(ratioDelta - phaseGps);

    metrics.chi_gps = chiGps;
    metrics.phase_gps = phaseGps;
    metrics.min_E_GPS = min(errorsGps);
    metrics.mean_E_GPS = mean(errorsGps);
    metrics.max_E_GPS = max(errorsGps);
end
