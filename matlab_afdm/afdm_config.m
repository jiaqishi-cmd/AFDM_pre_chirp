function config = afdm_config()
%AFDM_CONFIG 返回 AFDM 参数配置结构体，按功能分组。
%   config.modulation: 符号调制参数组
%   config.waveform: 波形参数组
%   config.channel: 信道参数组
%   config.simulation: 仿真参数组（循环配置）

    config.waveform.NumSubcarriers = 64;      % 子载波数量

    % 符号调制参数组
    config.modulation.M_mod = 16;             % 调制阶数（默认 16QAM）
    config.modulation.modType = 'qam';        % 调制类型 ('qam' 或 'psk')
    
    % 信道参数组
    config.channel.snr_db = 30;               % 信噪比（dB）
    config.channel.add_noise = true;         % 是否加噪 (true/false)
    config.channel.multipath = true;        % 是否多径 (true/false)
    taps = 3; max_delay = 3; k_max = 2.5;
    config.channel.delay_taps = randi([0, max_delay], 1, taps);          % 延迟采样数数组
    Doppler_taps = (2*rand(1,taps)-1) * k_max;  % 多普勒频移数组，范围 [-k_max, k_max]
    config.channel.doppler_freq = Doppler_taps / config.waveform.NumSubcarriers;        % 多普勒频移数组
    config.channel.chan_coef = (randn(1,taps) + 1i*randn(1,taps))/sqrt(2);           % 信道系数数组
    config.channel.chan_coef = config.channel.chan_coef / norm(config.channel.chan_coef);

    % 波形参数组
    config.waveform.CPPLength = 3;            % Chirp 循环前缀长度
    k_v = 1;
    config.waveform.c1 = (2*(floor(max(abs(Doppler_taps)))+k_v)+1) / (2*config.waveform.NumSubcarriers);                   % Post-chirp 参数
    config.waveform.c2 = sqrt(2) / (10 * config.waveform.NumSubcarriers);                    % Pre-chirp 参数
    
    % 仿真参数组（循环配置）
    config.simulation.enable_snr_loop = true;           % 是否启用SNR扫描循环
    config.simulation.snr_range = 0:5:30;               % SNR扫描范围（dB）
    config.simulation.num_trials = 1;                   % 每个SNR点的蒙特卡洛试验次数
    config.simulation.enable_montecarlo = false;        % 是否启用蒙特卡洛循环（单SNR点多次试验）

end
