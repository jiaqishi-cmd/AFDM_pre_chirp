function s = ifft_oversampled(xPre, os)
%IFFT_OVERSAMPLED Generate an oversampled time-domain waveform.

    M = numel(xPre);
    if os == 1
        s = ifft(xPre) * sqrt(M);
        return;
    end

    Mos = os * M;
    half = M / 2;
    Xos = zeros(Mos, 1);
    Xos(1:half) = xPre(1:half);
    Xos(end-half+1:end) = xPre(half+1:end);
    s = ifft(Xos) * sqrt(Mos);
end
