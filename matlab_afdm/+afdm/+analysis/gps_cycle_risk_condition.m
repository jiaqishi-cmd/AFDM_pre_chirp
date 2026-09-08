function result = gps_cycle_risk_condition( ...
        N, c1, relativeDelay, relativeDoppler, rotationalOrder)
%GPS_CYCLE_RISK_CONDITION Closed-form rational-GPS cycle condition.
%   Applies to integer delay/Doppler and the rational GPS phase magnitude
%   pi/2. rotationalOrder is 2 for BPSK and 4 for QPSK/square QAM.

    validateattributes(N, {'numeric'}, {'scalar', 'integer', 'positive'});
    validateattributes(relativeDelay, {'numeric'}, {'scalar', 'integer'});
    validateattributes(relativeDoppler, {'numeric'}, {'scalar', 'integer'});
    validateattributes(rotationalOrder, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});

    rawShift = relativeDoppler + 2 * N * c1 * relativeDelay;
    roundedShift = round(rawShift);
    if abs(rawShift - roundedShift) > 1e-9
        error('The closed-form condition requires an integer effective shift.');
    end

    eta = mod(-roundedShift, N);
    cycleCount = gcd(N, eta);
    cycleLength = N / cycleCount;
    phaseClosureResidue = mod( ...
        rotationalOrder * relativeDelay * eta, N);

    result.eta = eta;
    result.cycle_count = cycleCount;
    result.cycle_length = cycleLength;
    result.rotational_order = rotationalOrder;
    result.phase_closure_residue = phaseClosureResidue;
    result.has_cycle_excluding_zero = cycleCount > 1;
    result.predicts_compatible_cycle = ...
        result.has_cycle_excluding_zero && phaseClosureResidue == 0;
end
