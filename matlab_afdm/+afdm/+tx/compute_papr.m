function papr = compute_papr(signal)
%COMPUTE_PAPR Compute peak-to-average power ratio in dB.

    power = abs(signal).^2;
    papr = 10 * log10(max(power) / mean(power));
end
