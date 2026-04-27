function x_est = mmse_equalize(y_daft, H_eff, noise_var)
%MMSE_EQUALIZE 基于等效信道矩阵 H_eff 对 DAFT 符号进行 MMSE 均衡。
%   x_est = mmse_equalize(y_daft, H_eff, noise_var)

    [N, M] = size(H_eff);
    if N ~= M
        error('等效信道矩阵 H_eff 必须为方阵。');
    end
    if numel(y_daft) ~= N
        error('y_daft 长度必须等于 H_eff 的维度。');
    end

    W_mmse = (H_eff' * H_eff + noise_var * eye(N)) \ H_eff';
    x_est = W_mmse * y_daft(:);
end