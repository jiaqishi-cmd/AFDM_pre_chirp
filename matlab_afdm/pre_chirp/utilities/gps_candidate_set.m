function candidate_set = gps_candidate_set(numSubcarriers, numCandidates)
%GPS_CANDIDATE_SET Generate the GPS c2 candidate set Omega.
%   The W=2 case uses the original pair +/- 1/(4*m^2), with a small
%   constant value for m=0 to avoid the singular endpoint.

    validateattributes(numCandidates, {'numeric'}, {'scalar', 'integer', '>=', 2});

    candidate_set = zeros(numSubcarriers, numCandidates);

    for m_idx = 0:numSubcarriers-1
        if m_idx == 0
            base_val = 1 / (numSubcarriers^2);
        else
            base_val = 1 / (4 * m_idx^2);
        end

        for w = 1:numCandidates
            sign_val = (-1)^(w - 1);
            scale_val = ceil(w / 2);
            candidate_set(m_idx + 1, w) = sign_val * scale_val * base_val;
        end
    end
end
