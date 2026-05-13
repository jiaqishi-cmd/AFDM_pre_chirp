function xEst = mmse_equalize(yDaft, hEff, noiseVar)
%MMSE_EQUALIZE Equalize DAFT-domain symbols using the MMSE solution.

    [N, M] = size(hEff);
    if N ~= M
        error('H_eff must be square.');
    end
    if numel(yDaft) ~= N
        error('y_daft length must equal H_eff dimensions.');
    end

    wMmse = (hEff' * hEff + noiseVar * eye(N)) \ hEff';
    xEst = wMmse * yDaft(:);
end
