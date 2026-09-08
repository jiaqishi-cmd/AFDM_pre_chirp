function metrics = pairwise_diversity_metrics(psi, tolerance)
%PAIRWISE_DIVERSITY_METRICS Evaluate rank and coding-gain indicators.

    singularValues = svd(psi, 'econ');
    if nargin < 2 || isempty(tolerance)
        tolerance = max(size(psi)) * eps(max(singularValues));
    end

    numericalRank = sum(singularValues > tolerance);
    gramEigenvalues = sort(real(eig(psi' * psi)), 'descend');
    gramDeterminant = max(real(det(psi' * psi)), 0);
    positiveEigenvalues = gramEigenvalues(gramEigenvalues > tolerance^2);

    if isempty(singularValues)
        minSingularValue = 0;
        conditionNumber = Inf;
    else
        minSingularValue = singularValues(end);
        conditionNumber = singularValues(1) / max(minSingularValue, eps);
    end

    if isempty(positiveEigenvalues)
        pseudoDeterminant = 0;
    else
        pseudoDeterminant = prod(positiveEigenvalues);
    end

    metrics.rank = numericalRank;
    metrics.diversity_order = numericalRank;
    metrics.singular_values = singularValues;
    metrics.min_singular_value = minSingularValue;
    metrics.gram_eigenvalues = gramEigenvalues;
    metrics.gram_determinant = gramDeterminant;
    metrics.gram_pseudodeterminant = pseudoDeterminant;
    metrics.condition_number = conditionNumber;
    metrics.tolerance = tolerance;
end
