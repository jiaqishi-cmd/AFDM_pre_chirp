function xEst = equalize_symbols(yDaft, hEff, noiseVar, config)
%EQUALIZE_SYMBOLS Equalize DAFT-domain symbols with a selectable detector.

    detector = 'full_mmse';
    if nargin >= 4 && isstruct(config) && isfield(config, 'receiver') && ...
            isfield(config.receiver, 'detector') && ~isempty(config.receiver.detector)
        detector = config.receiver.detector;
    end

    switch lower(detector)
        case 'full_mmse'
            xEst = afdm.rx.mmse_equalize(yDaft, hEff, noiseVar);
        case 'diagonal_mmse'
            xEst = diagonal_mmse_equalize(yDaft, hEff, noiseVar);
        case 'banded_mmse'
            bandwidth = 1;
            if isfield(config.receiver, 'bandwidth') && ~isempty(config.receiver.bandwidth)
                bandwidth = config.receiver.bandwidth;
            end
            hSparse = apply_circular_band_mask(hEff, bandwidth);
            xEst = afdm.rx.mmse_equalize(yDaft, hSparse, noiseVar);
        case 'topk_mmse'
            topK = 4;
            if isfield(config.receiver, 'top_k') && ~isempty(config.receiver.top_k)
                topK = config.receiver.top_k;
            end
            hSparse = apply_row_topk_mask(hEff, topK);
            xEst = afdm.rx.mmse_equalize(yDaft, hSparse, noiseVar);
        case 'column_topk_mmse'
            topK = 4;
            if isfield(config.receiver, 'top_k') && ~isempty(config.receiver.top_k)
                topK = config.receiver.top_k;
            end
            hSparse = apply_column_topk_mask(hEff, topK);
            xEst = afdm.rx.mmse_equalize(yDaft, hSparse, noiseVar);
        otherwise
            error('Unsupported receiver detector: %s', detector);
    end
end

function xEst = diagonal_mmse_equalize(yDaft, hEff, noiseVar)
    h = diag(hEff);
    yDaft = yDaft(:);

    if numel(yDaft) ~= numel(h)
        error('y_daft length must equal H_eff dimensions.');
    end

    denom = abs(h).^2 + noiseVar;
    denom(denom == 0) = eps;
    xEst = conj(h) .* yDaft ./ denom;
end

function hSparse = apply_circular_band_mask(hEff, bandwidth)
    validateattributes(bandwidth, {'numeric'}, {'scalar', 'integer', 'nonnegative'});

    [N, M] = size(hEff);
    if N ~= M
        error('H_eff must be square.');
    end

    hSparse = zeros(size(hEff));
    for row = 1:N
        for offset = -bandwidth:bandwidth
            col = mod(row - 1 + offset, N) + 1;
            hSparse(row, col) = hEff(row, col);
        end
    end
end

function hSparse = apply_row_topk_mask(hEff, topK)
    validateattributes(topK, {'numeric'}, {'scalar', 'integer', 'positive'});

    [N, M] = size(hEff);
    topK = min(topK, M);
    hSparse = zeros(size(hEff));

    for row = 1:N
        [~, order] = sort(abs(hEff(row, :)).^2, 'descend');
        keep = order(1:topK);
        hSparse(row, keep) = hEff(row, keep);
    end
end

function hSparse = apply_column_topk_mask(hEff, topK)
    validateattributes(topK, {'numeric'}, {'scalar', 'integer', 'positive'});

    [N, M] = size(hEff);
    topK = min(topK, N);
    hSparse = zeros(size(hEff));

    for col = 1:M
        [~, order] = sort(abs(hEff(:, col)).^2, 'descend');
        keep = order(1:topK);
        hSparse(keep, col) = hEff(keep, col);
    end
end
