function [candidateSet, metadata] = gps_candidate_set_variant( ...
        numSubcarriers, variant, options)
%GPS_CANDIDATE_SET_VARIANT Build rational, irrational, or quantized GPS.

    if nargin < 3
        options = struct();
    end
    if nargin < 2 || isempty(variant)
        variant = 'rational';
    end

    precisionK = get_option(options, 'precision_k', 2);
    phaseBits = get_option(options, 'phase_bits', 8);
    irrationalScale = pi * 10^precisionK / floor(pi * 10^precisionK);

    switch lower(variant)
        case 'rational'
            phaseMagnitude = pi / 2;
        case 'irrational'
            phaseMagnitude = (pi / 2) * irrationalScale;
        case 'phase_quantized'
            phaseStep = 2 * pi / 2^phaseBits;
            phaseMagnitude = round(((pi / 2) * irrationalScale) / phaseStep) ...
                * phaseStep;
        case 'phase_offset'
            phaseMagnitude = pi / 2 + get_option(options, 'phase_offset', 0);
        otherwise
            error('Unsupported GPS variant: %s', variant);
    end

    candidateSet = zeros(numSubcarriers, 2);
    for m = 1:numSubcarriers-1
        c2Magnitude = phaseMagnitude / (2 * pi * m^2);
        candidateSet(m + 1, :) = [c2Magnitude, -c2Magnitude];
    end

    metadata.variant = variant;
    metadata.precision_k = precisionK;
    metadata.phase_bits = phaseBits;
    metadata.irrational_scale = irrationalScale;
    metadata.phase_magnitude = phaseMagnitude;
    metadata.phase_error_from_pi_over_2 = phaseMagnitude - pi / 2;
end

function value = get_option(options, name, defaultValue)
    if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
        value = options.(name);
    else
        value = defaultValue;
    end
end
