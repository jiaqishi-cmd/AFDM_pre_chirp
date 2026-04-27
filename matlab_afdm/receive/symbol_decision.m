function x_dec = symbol_decision(x_est, M_mod, modType)
%SYMBOL_DECISION 根据调制类型对均衡后的复数符号执行判决。
%   x_dec = symbol_decision(x_est, M_mod, modType)
%   modType: 'qam' 或 'psk'

    if nargin < 2 || isempty(M_mod)
        error('请输入有效的调制阶数 M_mod。');
    end
    if nargin < 3 || isempty(modType)
        modType = 'qam';
    end

    switch lower(modType)
        case 'qam'
            x_dec = qamdemod(x_est, M_mod, 'OutputType', 'bit', 'UnitAveragePower', true);
        case 'psk'
            x_dec = pskdemod(x_est, M_mod, 'OutputType', 'bit', 'PhaseOffset', pi/M_mod);
        otherwise
            error('不支持的调制类型: %s. 仅支持 ''qam'' 或 ''psk''.', modType);
    end
end