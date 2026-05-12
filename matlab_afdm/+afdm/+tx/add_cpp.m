function signalCpp = add_cpp(signal, cppLength, c1)
%ADD_CPP Add chirp cyclic prefix.

    if nargin < 3 || isempty(c1)
        c1 = 0;
    end

    if cppLength < 0 || cppLength > numel(signal)
        error('Invalid Chirp cyclic prefix length.');
    end

    N = numel(signal);
    signal = signal(:);

    if c1 == 0
        signalCpp = [signal(end-cppLength+1:end); signal];
    else
        signalCpp = zeros(N + cppLength, 1);
        signalCpp(cppLength+1:end) = signal;

        sTail = signal(N-cppLength+1:end);
        idxNeg = (-cppLength:-1).';
        cppPhase = exp(-1i * 2 * pi * c1 * (N^2 + 2 * N * idxNeg));
        signalCpp(1:cppLength) = sTail .* cppPhase;
    end
end
