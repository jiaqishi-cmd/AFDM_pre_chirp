function [delta_cont, delta_qpsk, delta_bpsk] = build_recursive_delta_from_gps(N, L, d_m, path_phase)
%BUILD_RECURSIVE_DELTA_FROM_GPS 根据 GPS 相位递推构造 delta。
%   目标是让 ratio_delta(z)=delta_z^2/(delta_{z-L}delta_{z+L})
%   尽量匹配 path_phase * chi_GPS(z,L)。这里使用循环松弛迭代求相位。

    d_m = d_m(:);
    z = (0:N-1).';
    z_minus = mod(z - L, N) + 1;
    z_plus = mod(z + L, N) + 1;
    target = path_phase * (d_m.^2 ./ (d_m(z_minus) .* d_m(z_plus)));
    target_phase = angle(target);

    theta = zeros(N, 1);
    maxIter = 200;
    for iter = 1:maxIter
        oldTheta = theta;
        for idx = 1:N
            neighbors = angle(exp(1i * (theta(z_minus(idx)) + theta(z_plus(idx)) + target_phase(idx))));
            theta(idx) = 0.5 * neighbors;
        end
        if norm(angle(exp(1i * (theta - oldTheta)))) < 1e-10
            break;
        end
    end

    delta_cont = exp(1i * theta);
    delta_qpsk = quantize_phase(delta_cont, [1, 1i, -1, -1i]);
    delta_bpsk = quantize_phase(delta_cont, [1, -1]);
end

function quantized = quantize_phase(values, alphabet)
    alphabet = alphabet(:).';
    quantized = zeros(size(values));
    for idx = 1:numel(values)
        [~, best] = min(abs(values(idx) - alphabet));
        quantized(idx) = alphabet(best);
    end
end
