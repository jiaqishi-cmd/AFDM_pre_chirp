function symbols = daft_demod(signal, numSubcarriers, c1, c2)
%DAFT_DEMOD Map a time-domain AFDM signal back to the DAFT domain.

    if nargin < 4 || isempty(c2)
        c2 = 0;
    end
    if nargin < 3 || isempty(c1)
        c1 = 0;
    end

    validateattributes(numSubcarriers, {'numeric'}, {'scalar', 'integer', 'positive'});
    N = numSubcarriers;
    signal = signal(:);
    if numel(signal) ~= N
        error('signal length must equal numSubcarriers.');
    end

    if isscalar(c2)
        c2Coeff = c2;
    else
        if numel(c2) ~= N
            error('c2 length must equal numSubcarriers.');
        end
        c2Coeff = c2(:);
    end

    n = (0:N-1).';
    l1Phase = exp(-1i * 2 * pi * c1 * (n.^2));
    xPre = signal .* l1Phase;

    xIfft = fft(xPre) / sqrt(N);

    l2Phase = exp(-1i * 2 * pi * c2Coeff .* (n.^2));
    symbols = xIfft .* l2Phase;
end
