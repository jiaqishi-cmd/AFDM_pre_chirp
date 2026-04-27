function signal_cpp = add_cpp(signal, cppLength, c1)
%ADD_CPP 添加 Chirp 循环前缀（CPP），可选应用 Chirp 相位。
%   signal_cpp = add_cpp(signal, cppLength, c1)
%   signal: 输入信号。
%   cppLength: Chirp 循环前缀长度。
%   c1: Post-chirp 参数（默认 0，为 0 时退化为普通 CPP）。
%
%   signal_cpp: 包含 Chirp 循环前缀的信号。

    if nargin < 3 || isempty(c1)
        c1 = 0;
    end
    
    if cppLength < 0 || cppLength > numel(signal)
        error('Invalid Chirp cyclic prefix length.');
    end
    
    N = numel(signal);
    signal = signal(:);
    
    if c1 == 0
        % 标准 CPP：直接复制信号尾部
        signal_cpp = [signal(end-cppLength+1:end); signal];
    else
        % 带 Chirp 相位的 CPP
        N_total = N + cppLength;
        signal_cpp = zeros(N_total, 1);
        signal_cpp(cppLength+1:end) = signal;
        
        % 提取信号尾部作为 CPP 基础
        s_tail = signal(N-cppLength+1:end);
        
        % 计算 Chirp 相位：针对 CPP 部分（负索引）
        idx_neg = (-cppLength:-1).';
        cpp_phase = exp(-1i * 2 * pi * c1 * (N^2 + 2*N*idx_neg));
        
        % 将相位应用到 Chirp 循环前缀
        signal_cpp(1:cppLength) = s_tail .* cpp_phase;
    end
end
