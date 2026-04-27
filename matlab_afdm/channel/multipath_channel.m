function y = multipath_channel(signal, config)
%MULTIPATH_CHANNEL 物理时域多径信道，支持延迟和多普勒效应。
%   y = multipath_channel(signal, config)
%   signal: 输入基带信号。
%   config: 包含 config.channel 的配置结构体。
%           config.channel.multipath: 是否启用多径。
%           config.channel.delay_taps: 延迟采样数数组（默认 [0]）。
%           config.channel.doppler_freq: 多普勒频移数组（默认 [0]）。
%           config.channel.chan_coef: 信道系数数组（默认 [1]）。
%           config.waveform.CPPLength: Chirp 循环前缀长度。
%
%   y: 通过信道的信号。

    if ~isfield(config, 'channel') || ~isfield(config.channel, 'multipath')
        error('config.channel.multipath not found.');
    end
    
    if ~config.channel.multipath
        % 多径不启用，直接返回原信号
        y = signal;
        return;
    end
    
    % 从配置中提取参数
    delay_taps = config.channel.delay_taps;
    doppler_freq = config.channel.doppler_freq;
    chan_coef = config.channel.chan_coef;
    cpp_len = config.waveform.CPPLength;
    
    % 验证参数长度一致性
    taps = numel(delay_taps);
    if numel(doppler_freq) ~= taps || numel(chan_coef) ~= taps
        error('delay_taps, doppler_freq, chan_coef 长度必须相同。');
    end
    
    signal = signal(:);
    N_total = numel(signal);
    
    % 时域卷积，考虑延迟和多普勒
    y = zeros(N_total, 1);
    for n_abs = 1:N_total
        t = n_abs - (cpp_len + 1);
        for p = 1:taps
            l = delay_taps(p);
            f = doppler_freq(p);
            if (n_abs - l) > 0
                val_s = signal(n_abs - l);
                h_val = chan_coef(p) * exp(-1i * 2 * pi * f * t);
                y(n_abs) = y(n_abs) + h_val * val_s;
            end
        end
    end
end
