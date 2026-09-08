function metrics = linear_detector_floor_metrics(hEff, MMod, modType, numFrames, seed)
%LINEAR_DETECTOR_FLOOR_METRICS Diagnose the zero-noise linear-detector floor.
%   Uses pinv(H)*H, the zero-noise limit of a full linear equalizer on a
%   rank-deficient channel, followed by the configured hard slicer.

    if nargin < 5 || isempty(seed)
        seed = 20260907;
    end
    if nargin < 4 || isempty(numFrames)
        numFrames = 10000;
    end

    [numRows, numCols] = size(hEff);
    if numRows ~= numCols
        error('hEff must be square.');
    end

    projector = pinv(hEff) * hEff;
    desiredPower = abs(diag(projector)).^2;
    interferencePower = sum(abs(projector).^2, 2) - desiredPower;
    rowSir = desiredPower ./ max(interferencePower, eps);

    oldRngState = rng;
    cleanup = onCleanup(@() rng(oldRngState)); %#ok<NASGU>
    rng(seed, 'twister');

    totalErrors = 0;
    totalBits = 0;
    collectMargins = MMod == 2 && strcmpi(modType, 'psk');
    if collectMargins
        decisionMargins = zeros(numCols * numFrames, 1);
    else
        decisionMargins = [];
    end
    for frameIdx = 1:numFrames
        [symbols, bits] = afdm.tx.random_data(numCols, MMod, modType);
        estimate = projector * symbols;
        if collectMargins
            firstIdx = (frameIdx - 1) * numCols + 1;
            lastIdx = frameIdx * numCols;
            decisionMargins(firstIdx:lastIdx) = real(conj(symbols) .* estimate);
        end
        decisions = afdm.rx.symbol_decision(estimate, MMod, modType);
        [errors, bitsInFrame] = afdm.rx.compute_bit_errors(decisions, bits);
        totalErrors = totalErrors + errors;
        totalBits = totalBits + bitsInFrame;
    end

    singularValues = svd(hEff);
    metrics.rank_h = rank(hEff, 1e-10);
    metrics.min_singular_h = singularValues(end);
    metrics.projector = projector;
    metrics.projector_error_fro = norm(projector - eye(numCols), 'fro');
    metrics.min_row_sir_db = 10 * log10(min(rowSir));
    metrics.median_row_sir_db = 10 * log10(median(rowSir));
    metrics.noiseless_bit_errors = totalErrors;
    metrics.total_bits = totalBits;
    metrics.noiseless_ber = totalErrors / totalBits;
    if collectMargins
        metrics.min_decision_margin = min(decisionMargins);
        metrics.p001_decision_margin = percentile_by_sort(decisionMargins, 0.001);
        metrics.p01_decision_margin = percentile_by_sort(decisionMargins, 0.01);
        metrics.nonpositive_margin_fraction = mean(decisionMargins <= 0);
    else
        metrics.min_decision_margin = NaN;
        metrics.p001_decision_margin = NaN;
        metrics.p01_decision_margin = NaN;
        metrics.nonpositive_margin_fraction = NaN;
    end
end

function value = percentile_by_sort(samples, probability)
    samples = sort(samples(:));
    index = max(1, min(numel(samples), ceil(probability * numel(samples))));
    value = samples(index);
end
