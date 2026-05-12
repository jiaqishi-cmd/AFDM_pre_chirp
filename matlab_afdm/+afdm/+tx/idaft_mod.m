function signal = idaft_mod(symbols, numSubcarriers, c1, c2)
%IDAFT_MOD Apply AFDM modulation. Supports scalar or vector c2.

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
    if isscalar(c2)
        c2Coeff = c2;
    else
        if numel(c2) ~= N
            error('c2 vector length must match numSubcarriers.');
        end
        c2Coeff = c2(:);
    end

    F = dftmtx(N);
    F = F ./ norm(F);

    L1 = diag(exp(-1i * 2 * pi * c1 * (n.^2)));
    L2 = diag(exp(-1i * 2 * pi * c2Coeff .* (n.^2)));

    A = L2 * F * L1;
    signal = A' * x;
end
