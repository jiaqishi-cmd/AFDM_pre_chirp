% 循环仿真脚本：支持SNR扫描和蒙特卡洛试验
% 此脚本演示 AFDM 系统的性能评估

% 添加子文件夹路径
addpath(genpath('channel'));
addpath(genpath('receive'));
addpath(genpath('transmitter'));

% 1. 获取配置
config = afdm_config();

% 2. 根据配置确定仿真模式
if config.simulation.enable_snr_loop
    % SNR 扫描循环模式
    snr_values = config.simulation.snr_range;
    num_trials = config.simulation.num_trials;
    
    fprintf('========== SNR 扫描循环仿真 ==========\n');
    fprintf('SNR 范围: %s dB\n', mat2str(snr_values));
    fprintf('每个 SNR 点的试验次数: %d\n\n', num_trials);
    
    % 初始化存储结果
    ber_results = zeros(length(snr_values), 1);
    papr_results = zeros(length(snr_values), 1);
    
    % SNR 循环
    for snr_idx = 1:length(snr_values)
        current_snr = snr_values(snr_idx);
        config.channel.snr_db = current_snr;
        
        total_err_bits = 0;
        total_bits = 0;
        total_papr = 0;
        
        % 蒙特卡洛试验循环
        for trial = 1:num_trials
            % 发射端处理
            [signal_cpp, papr, tx_bits] = afdm_tx_engine(config);
            total_papr = total_papr + papr;
            
            % 信道处理
            r_signal = multipath_channel(signal_cpp, config);
            r_signal = add_awgn(r_signal, config);
            
            % 接收端处理
            [~, err_bits, num_bits] = afdm_rx_engine(r_signal, config, tx_bits);
            
            % 累积错误统计
            total_err_bits = total_err_bits + err_bits;
            total_bits = total_bits + num_bits;
        end
        
        % 计算平均结果
        ber_results(snr_idx) = total_err_bits / total_bits;
        papr_results(snr_idx) = total_papr / num_trials;
        
        fprintf('SNR = %6.1f dB | BER = %.2e | 平均PAPR = %.2f dB\n', ...
            current_snr, ber_results(snr_idx), papr_results(snr_idx));
    end
    
    % 绘制结果
    figure('Name', 'SNR扫描结果');
    
    subplot(1, 2, 1);
    semilogy(snr_values, ber_results, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
    grid on;
    xlabel('SNR (dB)');
    ylabel('误码率 (BER)');
    title('AFDM 系统性能');
    
    subplot(1, 2, 2);
    plot(snr_values, papr_results, 'rs-', 'LineWidth', 2, 'MarkerSize', 8);
    grid on;
    xlabel('SNR (dB)');
    ylabel('PAPR (dB)');
    title('PAPR vs SNR');
    
    fprintf('\n========== 仿真完成 ==========\n');
    
elseif config.simulation.enable_montecarlo
    % 蒙特卡洛循环模式（单 SNR 点）
    num_trials = config.simulation.num_trials;
    fprintf('========== 蒙特卡洛循环仿真 (SNR = %.1f dB) ==========\n', config.channel.snr_db);
    fprintf('试验次数: %d\n\n', num_trials);
    
    total_err_bits = 0;
    total_bits = 0;
    total_papr = 0;
    
    for trial = 1:num_trials
        % 发射端处理
        [signal_cpp, papr, tx_bits] = afdm_tx_engine(config);
        total_papr = total_papr + papr;
        
        % 信道处理
        r_signal = multipath_channel(signal_cpp, config);
        r_signal = add_awgn(r_signal, config);
        
        % 接收端处理
        [~, err_bits, num_bits] = afdm_rx_engine(r_signal, config, tx_bits);
        
        % 累积错误统计
        total_err_bits = total_err_bits + err_bits;
        total_bits = total_bits + num_bits;
        
        if mod(trial, 10) == 0 || trial == 1
            fprintf('试验 %4d/%d 完成\n', trial, num_trials);
        end
    end
    
    % 计算平均结果
    ber = total_err_bits / total_bits;
    avg_papr = total_papr / num_trials;
    
    fprintf('\n========== 结果统计 ==========\n');
    fprintf('总比特数: %d\n', total_bits);
    fprintf('总误比特数: %d\n', total_err_bits);
    fprintf('误码率 (BER): %.2e\n', ber);
    fprintf('平均 PAPR: %.2f dB\n', avg_papr);
    fprintf('========== 仿真完成 ==========\n');
    
else
    % 单次仿真模式
    fprintf('========== 单次仿真 ==========\n');
    
    % 发射端处理
    [signal_cpp, papr, tx_bits] = afdm_tx_engine(config);
    fprintf('发射信号 PAPR: %.2f dB\n', papr);
    
    % 信道处理
    r_signal = multipath_channel(signal_cpp, config);
    r_signal = add_awgn(r_signal, config);
    
    % 接收端处理
    [x_dec, err_bits, total_bits] = afdm_rx_engine(r_signal, config, tx_bits);
    
    % 结果显示
    ber = err_bits / total_bits;
    fprintf('总比特数: %d\n', total_bits);
    fprintf('误比特数: %d\n', err_bits);
    fprintf('误码率 (BER): %.2e\n', ber);
    fprintf('========== 仿真完成 ==========\n');
end
