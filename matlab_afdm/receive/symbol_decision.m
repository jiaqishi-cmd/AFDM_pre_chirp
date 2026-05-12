function x_dec = symbol_decision(x_est, M_mod, modType)
%SYMBOL_DECISION 根据调制类型对均衡后的复数符号执行判决。
%   x_dec = symbol_decision(x_est, M_mod, modType)
%   modType: 'qam' 或 'psk'

    if nargin < 2 || isempty(M_mod)
        error('M_mod must be provided.');
    end
    if nargin < 3 || isempty(modType)
        modType = 'qam';
    end

    switch lower(modType)
        case 'qam'
            x_dec = qamdemod(x_est, M_mod, 'OutputType', 'bit', 'UnitAveragePower', true);
        case 'psk'
            x_dec = pskdemod(x_est, M_mod, pi/M_mod, 'OutputType', 'bit');
        otherwise
            error('Unsupported modulation type: %s. Supported types are qam and psk.', modType);
    end
end
