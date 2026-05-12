function candidateSet = gps_candidate_set(numSubcarriers, numCandidates)
%GPS_CANDIDATE_SET Generate the GPS c2 candidate set Omega.
%   The W=2 case uses the original pair +/- 1/(4*m^2), with a small
%   constant value for m=0 to avoid the singular endpoint.

    validateattributes(numCandidates, {'numeric'}, {'scalar', 'integer', '>=', 2});

    candidateSet = zeros(numSubcarriers, numCandidates);
    for mIdx = 0:numSubcarriers-1
        if mIdx == 0
            baseVal = 1 / (numSubcarriers^2);
        else
            baseVal = 1 / (4 * mIdx^2);
        end

        for w = 1:numCandidates
            signVal = (-1)^(w - 1);
            scaleVal = ceil(w / 2);
            candidateSet(mIdx + 1, w) = signVal * scaleVal * baseVal;
        end
    end
end
