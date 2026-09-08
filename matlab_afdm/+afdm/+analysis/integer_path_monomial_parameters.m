function path = integer_path_monomial_parameters( ...
        N, c1, c2, delay, doppler)
%INTEGER_PATH_MONOMIAL_PARAMETERS Closed-form integer path representation.
%   H*e_q = coefficient(q)*e_permutation(q), using one-based MATLAB indices.

    validateattributes(N, {'numeric'}, {'scalar', 'integer', 'positive'});
    validateattributes(delay, {'numeric'}, {'scalar', 'integer', 'nonnegative'});
    validateattributes(doppler, {'numeric'}, {'scalar', 'integer'});

    if isscalar(c2)
        c2 = repmat(c2, N, 1);
    elseif numel(c2) ~= N
        error('c2 must be scalar or length N.');
    else
        c2 = c2(:);
    end

    effectiveOffset = doppler + 2 * N * c1 * delay;
    roundedOffset = round(effectiveOffset);
    if abs(effectiveOffset - roundedOffset) > 1e-9
        error('Integer-path monomial form requires integer 2*N*c1*delay+doppler.');
    end

    eta = mod(-roundedOffset, N);
    q = (0:N-1).';
    p = mod(q + eta, N);
    baseCoefficient = exp(1i * 2 * pi / N * ...
        (N * c1 * delay^2 - q * delay));
    preChirpPhase = exp(-1i * 2 * pi * c2 .* (q.^2));
    coefficient = baseCoefficient .* ...
        preChirpPhase(p + 1) .* conj(preChirpPhase(q + 1));

    path.eta = eta;
    path.input_index_zero_based = q;
    path.output_index_zero_based = p;
    path.permutation = p + 1;
    path.coefficient = coefficient;
    path.base_coefficient = baseCoefficient;
    path.pre_chirp_phase = preChirpPhase;
    path.delay = delay;
    path.doppler = doppler;
end
