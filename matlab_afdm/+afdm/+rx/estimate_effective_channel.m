function hEff = estimate_effective_channel(N, c1, c2, gains, delays, dopplerFreq)
%ESTIMATE_EFFECTIVE_CHANNEL Estimate the AFDM-domain effective channel matrix.

    P = length(gains);
    if length(delays) ~= P || length(dopplerFreq) ~= P
        error('gains, delays, and doppler_freq must have the same length.');
    end
    if isscalar(c2)
        c2 = repmat(c2, N, 1);
    elseif numel(c2) ~= N
        error('c2 must be scalar or a vector of length N.');
    else
        c2 = c2(:);
    end

    hTime = zeros(N, N);
    nVec = (0:N-1).';
    identityN = eye(N);

    for idx = 1:P
        hPath = gains(idx);
        delay = delays(idx);
        nu = dopplerFreq(idx);

        piMat = circshift(identityN, delay, 1);
        phaseDoppler = exp(-1i * 2 * pi * nu * nVec);
        deltaMat = diag(phaseDoppler);

        gammaDiag = ones(N, 1);
        for n = 0:N-1
            if n < delay
                phaseCpp = -1i * 2 * pi * c1 * (N^2 - 2 * N * (delay - n));
                gammaDiag(n+1) = exp(phaseCpp);
            end
        end
        gammaMat = diag(gammaDiag);

        hTime = hTime + hPath * gammaMat * deltaMat * piMat;
    end

    hEff = zeros(N, N);
    n = nVec;
    phasePostDemod = exp(-1i * 2 * pi * c1 * (n.^2));
    phasePreMod = exp(1i * 2 * pi * c2(:) .* (n.^2));
    phasePreDemod = exp(-1i * 2 * pi * c2(:) .* (n.^2));
    phasePostMod = exp(1i * 2 * pi * c1 * (n.^2));

    for k = 1:N
        delta = zeros(N, 1);
        delta(k) = 1;

        xPre = delta .* phasePreMod;
        xTime = ifft(xPre) * sqrt(N);
        xTime = xTime .* phasePostMod;

        yTime = hTime * xTime;

        yPost = yTime .* phasePostDemod;
        xIfft = fft(yPost) / sqrt(N);
        yDaft = xIfft .* phasePreDemod;

        hEff(:, k) = yDaft;
    end
end
