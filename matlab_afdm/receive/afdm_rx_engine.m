function [x_dec, err_symbols, total_symbols] = afdm_rx_engine(r_signal, config, x_ref)
%AFDM_RX_ENGINE 接收端处理：去 CPP、DAFT 解调、MMSE 均衡、判决和误符号统计。
%   [x_dec, err_symbols, total_symbols] = afdm_rx_engine(r_signal, config, x_ref)
%
%   参数均从 config 中读取，包括信噪比 snr_db。

    % 1. 移除 CPP
    r_data = afdm.rx.remove_cpp(r_signal, config.waveform.CPPLength);

    % 2. DAFT 解调
    y_daft = afdm.rx.daft_demod(r_data, config.waveform.NumSubcarriers, config.waveform.c1, config.waveform.c2);

    % 3. 估计当前等效信道矩阵
    H_eff = afdm.rx.estimate_effective_channel(config.waveform.NumSubcarriers, config.waveform.c1, config.waveform.c2, config.channel.chan_coef, config.channel.delay_taps, config.channel.doppler_freq);

    % 4. 根据 SNR 计算噪声方差
    noise_var = 1 / (10^(config.channel.snr_db / 10));

    % 5. MMSE 均衡
    x_est = afdm.rx.mmse_equalize(y_daft, H_eff, noise_var);

    % 6. 判决
    x_dec = afdm.rx.symbol_decision(x_est, config.modulation.M_mod, config.modulation.modType);

    % 7. 统计误符号
    [err_symbols, total_symbols] = afdm.rx.compute_bit_errors(x_dec, x_ref);
end
