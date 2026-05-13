function dmin = dmin_over_delta_set(Htot, deltaSet)
%DMIN_OVER_DELTA_SET Compute normalized minimum distance over finite deltas.

    Y = Htot * deltaSet.delta;
    metric = sum(abs(Y).^2, 1) ./ deltaSet.norm2;
    [bestMetric, bestIdx] = min(real(metric));

    dmin.dmin2 = bestMetric;
    dmin.best_index = bestIdx;
    dmin.best_delta = deltaSet.delta(:, bestIdx);
    dmin.best_delta_type = deltaSet.type(bestIdx);
    dmin.best_delta_index = deltaSet.index(bestIdx);
end
