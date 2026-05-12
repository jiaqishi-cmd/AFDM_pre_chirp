function dmin = compute_dmin_over_delta_set(Htot, delta_set)
%COMPUTE_DMIN_OVER_DELTA_SET 计算有限 delta 集合上的归一化最小距离。
%   metric(delta)=||Htot*delta||^2/||delta||^2。

    Y = Htot * delta_set.delta;
    metric = sum(abs(Y).^2, 1) ./ delta_set.norm2;
    [bestMetric, bestIdx] = min(real(metric));

    dmin.dmin2 = bestMetric;
    dmin.best_index = bestIdx;
    dmin.best_delta = delta_set.delta(:, bestIdx);
    dmin.best_delta_type = delta_set.type(bestIdx);
    dmin.best_delta_index = delta_set.index(bestIdx);
end
