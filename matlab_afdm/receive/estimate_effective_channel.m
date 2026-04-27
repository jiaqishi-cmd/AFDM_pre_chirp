function H_eff = estimate_effective_channel(N, c1, c2, gains, delays, doppler_freq)
%ESTIMATE_EFFECTIVE_CHANNEL 估计 AFDM 等效信道矩阵 H_eff。
%   H_eff = estimate_effective_channel(N, c1, c2, gains, delays, doppler_freq)
%   N: 子载波/符号数量
%   c1: AFDM 线性调频参数
%   c2: 标量或向量形式的 c2 参数（标量或长度 N）
%   gains: 路径增益向量
%   delays: 路径整数延迟向量
%   doppler_freq: 路径归一化多普勒频移向量

    P = length(gains);
    if length(delays) ~= P || length(doppler_freq) ~= P
        error('增益、延迟和多普勒向量的长度必须一致');
    end
    if isscalar(c2)
        c2 = repmat(c2, N, 1);
    elseif numel(c2) ~= N
        error('c2 必须是标量或长度为 N 的向量');
    else
        c2 = c2(:);
    end

    % 生成 AFDM 时域信道矩阵 H_time
    H_time = zeros(N, N);
    n_vec = (0:N-1).';
    I_N = eye(N);

    for i = 1:P
        h_i = gains(i);
        l_i = delays(i);
        nu_i = doppler_freq(i);

        Pi_mat = circshift(I_N, l_i, 1);
        phase_doppler = exp(-1i * 2 * pi * nu_i * n_vec);
        Delta_mat = diag(phase_doppler);

        gamma_diag = ones(N, 1);
        for n = 0:N-1
            if n < l_i
                phase_cpp = -1i * 2 * pi * c1 * (N^2 - 2 * N * (l_i - n));
                gamma_diag(n+1) = exp(phase_cpp);
            end
        end
        Gamma_mat = diag(gamma_diag);

        H_time = H_time + h_i * Gamma_mat * Delta_mat * Pi_mat;
    end

    % 逐列构造等效信道矩阵 H_eff
    H_eff = zeros(N, N);
    n = n_vec;
    phase_post_demod = exp(-1i * 2 * pi * c1 * (n.^2));
    phase_pre_mod = exp(1i * 2 * pi * c2(:) .* (n.^2));
    phase_pre_demod = exp(-1i * 2 * pi * c2(:) .* (n.^2));
    phase_post_mod = exp(1i * 2 * pi * c1 * (n.^2));

    for k = 1:N
        delta = zeros(N, 1);
        delta(k) = 1;

        x_pre = delta .* phase_pre_mod;
        x_time = ifft(x_pre) * sqrt(N);
        x_time = x_time .* phase_post_mod;

        y_time = H_time * x_time;

        y_post = y_time .* phase_post_demod;
        x_ifft = fft(y_post) / sqrt(N);
        y_daft = x_ifft .* phase_pre_demod;

        H_eff(:, k) = y_daft;
    end
end
