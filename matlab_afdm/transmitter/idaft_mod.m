function signal = idaft_mod(symbols, numSubcarriers, c1, c2)
%IDAFT_MOD 对符号进行 AFDM 调制。支持 c2 为标量或向量。
%   signal = idaft_mod(symbols, numSubcarriers, c1, c2)
%   symbols: 输入调制符号向量。
%   numSubcarriers: 子载波数量。
%   c1: AFDM Chirp 参数（Post-chirp，标量，默认 0）。
%   c2: AFDM Chirp 参数（Pre-chirp）。
%       - 标量：所有子载波使用相同 c2（默认 0）
%       - 向量：每个子载波使用不同的 c2(n)
%
%   signal: AFDM 调制后的基带信号。

    if nargin < 4 || isempty(c2)
        c2 = 0;
    end
    if nargin < 3 || isempty(c1)
        c1 = 0;
    end

    N = numSubcarriers;
    
    if numel(symbols) > N
        error('Number of symbols exceeds number of subcarriers.');
    end
    
    x = symbols(:);
    if numel(x) < N
        x = [x; zeros(N - numel(x), 1)];
    end
    
    n = (0:N-1).';
    
    % 统一处理 c2 为标量或向量
    if isscalar(c2)
        c2_coeff = c2;
    else
        if numel(c2) ~= N
            error('c2 vector length must match numSubcarriers.');
        end
        c2_coeff = c2(:);
    end
    
    % 矩阵形式：A = L2 * F * L1, signal = A' * x
    F = dftmtx(N);
    F = F ./ norm(F);
    
    L1 = diag(exp(-1i * 2 * pi * c1 * (n.^2)));
    L2 = diag(exp(-1i * 2 * pi * c2_coeff .* (n.^2)));
    
    A = L2 * F * L1;
    signal = A' * x;
end
