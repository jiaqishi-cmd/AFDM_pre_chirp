function xDec = symbol_decision(xEst, MMod, modType)
%SYMBOL_DECISION Make hard bit decisions for QAM or PSK symbols.

    if nargin < 2 || isempty(MMod)
        error('M_mod must be provided.');
    end
    if nargin < 3 || isempty(modType)
        modType = 'qam';
    end

    switch lower(modType)
        case 'qam'
            xDec = qamdemod(xEst, MMod, 'OutputType', 'bit', 'UnitAveragePower', true);
        case 'psk'
            xDec = pskdemod(xEst, MMod, pi/MMod, 'OutputType', 'bit');
        otherwise
            error('Unsupported modulation type: %s. Supported types are qam and psk.', modType);
    end
end
