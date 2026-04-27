function symbols = daft_demod(signal, numSubcarriers, c1, c2)
%DAFT_DEMOD AFDM 解调：时域到 DAFT 域映射。
%   symbols = daft_demod(signal, numSubcarriers, c1, c2)
%   signal: 时域接收信号。
%   numSubcarriers: 子载波数量。
%   c1: AFDM Chirp 参数（Post-chirp，标量）。
%   c2: AFDM Chirp 参数（Pre-chirp，标量或向量）。
%
%   symbols: DAFT 域符号向量。

    if nargin < 4 || isempty(c2)
        c2 = 0;
    end
    if nargin < 3 || isempty(c1)
        c1 = 0;
    end

    N = numSubcarriers;
    if numel(signal) < N
        signal = [signal(:); zeros(N - numel(signal), 1)];
    else
        signal = signal(:);
    end

    if isscalar(c2)
        c2_coeff = c2;
    else
        if numel(c2) ~= N
            error('c2 length must equal numSubcarriers.');
        end
        c2_coeff = c2(:);
    end

    n = (0:N-1).';
    
    L1_phase = exp(-1i * 2 * pi * c1 * (n.^2));
    x_pre = signal .* L1_phase;

    x_ifft = fft(x_pre) / sqrt(N);

    L2_phase = exp(-1i * 2 * pi * c2_coeff .* (n.^2));
    symbols = x_ifft .* L2_phase;
end
