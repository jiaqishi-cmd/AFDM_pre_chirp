function deltaSet = generate_sparse_difference_set(N, MMod, modType, maxWeight)
%GENERATE_SPARSE_DIFFERENCE_SET Enumerate low-weight constellation errors.

    if nargin < 4 || isempty(maxWeight)
        maxWeight = min(N, 2);
    end
    if nargin < 3 || isempty(modType)
        modType = 'psk';
    end

    constellation = build_constellation(MMod, modType);
    differenceAlphabet = unique(round_complex( ...
        constellation(:) - constellation(:).', 12));
    differenceAlphabet = differenceAlphabet(abs(differenceAlphabet) > 1e-12);

    numEvents = 0;
    for weight = 1:maxWeight
        numEvents = numEvents + nchoosek(N, weight) * numel(differenceAlphabet)^weight;
    end
    delta = zeros(N, numEvents);
    weightLabel = zeros(numEvents, 1);
    cursor = 0;

    for weight = 1:maxWeight
        positionSets = nchoosek(1:N, weight);
        numAssignments = numel(differenceAlphabet)^weight;
        for positionIdx = 1:size(positionSets, 1)
            positions = positionSets(positionIdx, :);
            for assignment = 0:numAssignments-1
                digits = base_digits(assignment, numel(differenceAlphabet), weight);
                cursor = cursor + 1;
                delta(positions, cursor) = differenceAlphabet(digits + 1);
                weightLabel(cursor) = weight;
            end
        end
    end

    deltaSet.delta = delta;
    deltaSet.weight = weightLabel;
    deltaSet.norm2 = sum(abs(delta).^2, 1);
    deltaSet.difference_alphabet = differenceAlphabet;
    deltaSet.M_mod = MMod;
    deltaSet.mod_type = modType;
end

function constellation = build_constellation(MMod, modType)
    switch lower(modType)
        case 'psk'
            constellation = pskmod((0:MMod-1).', MMod, pi / MMod);
        case 'qam'
            constellation = qammod((0:MMod-1).', MMod, 'UnitAveragePower', true);
        otherwise
            error('Unsupported modulation type: %s', modType);
    end
end

function values = round_complex(values, digits)
    scale = 10^digits;
    values = round(real(values) * scale) / scale + ...
        1i * round(imag(values) * scale) / scale;
end

function digits = base_digits(value, base, width)
    digits = zeros(1, width);
    for idx = 1:width
        digits(idx) = mod(value, base);
        value = floor(value / base);
    end
end
