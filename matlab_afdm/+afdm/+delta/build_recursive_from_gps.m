function [deltaCont, deltaQpsk, deltaBpsk] = build_recursive_from_gps(N, L, dM, pathPhase)
%BUILD_RECURSIVE_FROM_GPS Build recursive delta candidates from GPS phase.

    dM = dM(:);
    z = (0:N-1).';
    zMinus = mod(z - L, N) + 1;
    zPlus = mod(z + L, N) + 1;
    target = pathPhase * (dM.^2 ./ (dM(zMinus) .* dM(zPlus)));
    targetPhase = angle(target);

    theta = zeros(N, 1);
    maxIter = 200;
    for iter = 1:maxIter
        oldTheta = theta;
        for idx = 1:N
            neighbors = angle(exp(1i * (theta(zMinus(idx)) + theta(zPlus(idx)) + targetPhase(idx))));
            theta(idx) = 0.5 * neighbors;
        end
        if norm(angle(exp(1i * (theta - oldTheta)))) < 1e-10
            break;
        end
    end

    deltaCont = exp(1i * theta);
    deltaQpsk = quantize_phase(deltaCont, [1, 1i, -1, -1i]);
    deltaBpsk = quantize_phase(deltaCont, [1, -1]);
end

function quantized = quantize_phase(values, alphabet)
    alphabet = alphabet(:).';
    quantized = zeros(size(values));
    for idx = 1:numel(values)
        [~, best] = min(abs(values(idx) - alphabet));
        quantized(idx) = alphabet(best);
    end
end
