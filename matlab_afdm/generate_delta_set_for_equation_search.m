function delta_set = generate_delta_set_for_equation_search(N, V, numRandomBpsk, numRandomQpsk, seed)
%GENERATE_DELTA_SET_FOR_EQUATION_SEARCH Build structured and random deltas.

    if nargin < 3 || isempty(numRandomBpsk)
        numRandomBpsk = 1000;
    end
    if nargin < 4 || isempty(numRandomQpsk)
        numRandomQpsk = 1000;
    end
    if nargin < 5 || isempty(seed)
        seed = 20260427;
    end

    M = N / V;
    if abs(M - round(M)) > 0
        error('N must be divisible by V.');
    end

    delta_set = struct('name', {}, 'delta', {});
    m = (0:N-1).';

    delta_set(end+1) = make_delta('delta1_all_ones', ones(N, 1));
    delta_set(end+1) = make_delta('delta2_alternating', (-1) .^ m);
    delta_set(end+1) = make_delta('delta3_group_bpsk', group_values(N, M, [1, -1, 1, -1]));
    delta_set(end+1) = make_delta('delta4_halves_same', [ones(N/2, 1); ones(N/2, 1)]);
    delta_set(end+1) = make_delta('delta5_halves_opposite', [ones(N/2, 1); -ones(N/2, 1)]);
    delta_set(end+1) = make_delta('delta6_group_qpsk', group_values(N, M, [1, 1i, -1, -1i]));
    delta_set(end+1) = make_delta('delta7_group_qpsk_alt', group_values(N, M, [1, -1i, -1, 1i]));

    oldRngState = rng;
    cleanup = onCleanup(@() rng(oldRngState));
    rng(seed, 'twister');

    for idx = 1:numRandomBpsk
        values = 2 * randi([0, 1], N, 1) - 1;
        delta_set(end+1) = make_delta(sprintf('delta8_random_bpsk_%04d', idx), values); %#ok<AGROW>
    end

    qpsk = [1, -1, 1i, -1i].';
    for idx = 1:numRandomQpsk
        values = qpsk(randi(numel(qpsk), N, 1));
        delta_set(end+1) = make_delta(sprintf('delta9_random_qpsk_%04d', idx), values); %#ok<AGROW>
    end
end

function item = make_delta(name, delta)
    item.name = name;
    item.delta = delta(:);
    if any(item.delta == 0)
        error('Delta entries must be nonzero.');
    end
end

function values = group_values(N, M, groupPattern)
    values = zeros(N, 1);
    for group_id = 1:numel(groupPattern)
        values((group_id - 1) * M + 1:group_id * M) = groupPattern(group_id);
    end
end
