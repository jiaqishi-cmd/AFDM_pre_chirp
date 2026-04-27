function r_data = remove_cpp(r_signal, cpp_len)
%REMOVE_CPP 去除 AFDM 接收信号中的 CPP 前缀。
%   r_data = remove_cpp(r_signal, cpp_len)
%   r_signal: 接收时域信号。
%   cpp_len: CPP 长度。

    if nargin < 2 || isempty(cpp_len)
        cpp_len = 0;
    end
    if cpp_len < 0 || cpp_len > numel(r_signal)
        error('CPP 长度必须在 0 到信号长度之间。');
    end

    r_data = r_signal(cpp_len + 1 : end);
end