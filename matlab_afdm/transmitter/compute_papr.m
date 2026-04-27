function papr = compute_papr(signal)
%COMPUTE_PAPR 计算 PAPR 值（峰均功率比）。
    power = abs(signal).^2;
    papr = 10 * log10(max(power) / mean(power));
end
