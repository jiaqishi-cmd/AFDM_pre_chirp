function x_est = equalize_symbols(y_daft, H_eff, noise_var, config)
%EQUALIZE_SYMBOLS Equalize DAFT-domain symbols with a selectable detector.

    detector = 'full_mmse';
    if nargin >= 4 && isstruct(config) && isfield(config, 'receiver') && ...
            isfield(config.receiver, 'detector') && ~isempty(config.receiver.detector)
        detector = config.receiver.detector;
    end

    switch lower(detector)
        case 'full_mmse'
            x_est = mmse_equalize(y_daft, H_eff, noise_var);
        case 'diagonal_mmse'
            x_est = diagonal_mmse_equalize(y_daft, H_eff, noise_var);
        case 'banded_mmse'
            bandwidth = 1;
            if isfield(config.receiver, 'bandwidth') && ~isempty(config.receiver.bandwidth)
                bandwidth = config.receiver.bandwidth;
            end
            H_sparse = apply_circular_band_mask(H_eff, bandwidth);
            x_est = mmse_equalize(y_daft, H_sparse, noise_var);
        case 'topk_mmse'
            topK = 4;
            if isfield(config.receiver, 'top_k') && ~isempty(config.receiver.top_k)
                topK = config.receiver.top_k;
            end
            H_sparse = apply_row_topk_mask(H_eff, topK);
            x_est = mmse_equalize(y_daft, H_sparse, noise_var);
        case 'column_topk_mmse'
            topK = 4;
            if isfield(config.receiver, 'top_k') && ~isempty(config.receiver.top_k)
                topK = config.receiver.top_k;
            end
            H_sparse = apply_column_topk_mask(H_eff, topK);
            x_est = mmse_equalize(y_daft, H_sparse, noise_var);
        otherwise
            error('Unsupported receiver detector: %s', detector);
    end
end

function x_est = diagonal_mmse_equalize(y_daft, H_eff, noise_var)
    h = diag(H_eff);
    y_daft = y_daft(:);

    if numel(y_daft) ~= numel(h)
        error('y_daft length must equal H_eff dimensions.');
    end

    denom = abs(h).^2 + noise_var;
    denom(denom == 0) = eps;
    x_est = conj(h) .* y_daft ./ denom;
end

function H_sparse = apply_circular_band_mask(H_eff, bandwidth)
    validateattributes(bandwidth, {'numeric'}, {'scalar', 'integer', 'nonnegative'});

    [N, M] = size(H_eff);
    if N ~= M
        error('H_eff must be square.');
    end

    H_sparse = zeros(size(H_eff));
    for row = 1:N
        for offset = -bandwidth:bandwidth
            col = mod(row - 1 + offset, N) + 1;
            H_sparse(row, col) = H_eff(row, col);
        end
    end
end

function H_sparse = apply_row_topk_mask(H_eff, topK)
    validateattributes(topK, {'numeric'}, {'scalar', 'integer', 'positive'});

    [N, M] = size(H_eff);
    topK = min(topK, M);
    H_sparse = zeros(size(H_eff));

    for row = 1:N
        [~, order] = sort(abs(H_eff(row, :)).^2, 'descend');
        keep = order(1:topK);
        H_sparse(row, keep) = H_eff(row, keep);
    end
end

function H_sparse = apply_column_topk_mask(H_eff, topK)
    validateattributes(topK, {'numeric'}, {'scalar', 'integer', 'positive'});

    [N, M] = size(H_eff);
    topK = min(topK, N);
    H_sparse = zeros(size(H_eff));

    for col = 1:M
        [~, order] = sort(abs(H_eff(:, col)).^2, 'descend');
        keep = order(1:topK);
        H_sparse(keep, col) = H_eff(keep, col);
    end
end
