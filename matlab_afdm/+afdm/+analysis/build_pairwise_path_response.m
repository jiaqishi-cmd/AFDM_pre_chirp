function psi = build_pairwise_path_response(pathMatrices, deltaX)
%BUILD_PAIRWISE_PATH_RESPONSE Build the path-response matrix for one error event.
%   The p-th column is H_p * deltaX. The rank of this matrix is the
%   pairwise diversity order under independent Rayleigh path gains.

    deltaX = deltaX(:);

    if iscell(pathMatrices)
        numPaths = numel(pathMatrices);
        if numPaths == 0
            error('pathMatrices must contain at least one path matrix.');
        end
        [numRows, numCols] = size(pathMatrices{1});
        psi = zeros(numRows, numPaths, 'like', pathMatrices{1});
        for pathIdx = 1:numPaths
            hPath = pathMatrices{pathIdx};
            if ~isequal(size(hPath), [numRows, numCols])
                error('All path matrices must have the same dimensions.');
            end
            if numel(deltaX) ~= numCols
                error('deltaX length must match the path-matrix column count.');
            end
            psi(:, pathIdx) = hPath * deltaX;
        end
        return;
    end

    if ndims(pathMatrices) ~= 3
        error('pathMatrices must be a cell array or an N-by-N-by-P array.');
    end
    [numRows, numCols, numPaths] = size(pathMatrices);
    if numel(deltaX) ~= numCols
        error('deltaX length must match the path-matrix column count.');
    end

    psi = zeros(numRows, numPaths, 'like', pathMatrices);
    for pathIdx = 1:numPaths
        psi(:, pathIdx) = pathMatrices(:, :, pathIdx) * deltaX;
    end
end
