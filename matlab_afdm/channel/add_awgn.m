function y = add_awgn(signal, config)
%ADD_AWGN 添加高斯白噪声。
%   y = add_awgn(signal, config)
%   signal: 输入信号。
%   config: 包含 config.channel 的配置结构体。
%           config.channel.snr_db: 信噪比（dB）。
%           config.channel.add_noise: 是否加噪（布尔值）。
%
%   y: 加噪后的信号。

    if ~isfield(config, 'channel') || ~isfield(config.channel, 'snr_db')
        error('config.channel.snr_db not found.');
    end
    
    signal = complex(signal);
    snr_db = config.channel.snr_db;
    
    if ~config.channel.add_noise
        % 不加噪，直接返回原信号
        y = signal;
        return;
    end
    
    power = mean(abs(signal).^2);
    snr = 10^(snr_db / 10);
    noise_power = power / snr;
    
    noise = sqrt(noise_power / 2) * (randn(size(signal)) + 1i * randn(size(signal)));
    y = signal + noise;

end
