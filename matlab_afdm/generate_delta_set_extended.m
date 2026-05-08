function delta_set = generate_delta_set_extended(N, V, L, d_gps, path_phase, options)
%GENERATE_DELTA_SET_EXTENDED 生成结构化、随机和递推 delta 集合。
%   delta 只用于 Bemani 关键等式和 Phi(delta) 搜索，不要求一定来自
%   真实 QAM 符号差分。所有元素均非零。

    if nargin < 6
        options = struct();
    end

    M = N / V;
    if abs(M - round(M)) > 0
        error('N must be divisible by V.');
    end

    maxExhaustiveBpsk = get_option(options, 'max_exhaustive_bpsk', 4096);
    numRandomBpsk = get_option(options, 'num_random_bpsk', 2000);
    numRandomQpsk = get_option(options, 'num_random_qpsk', 2000);
    seed = get_option(options, 'seed', 20260427);

    delta_set = struct('name', {}, 'delta', {});
    m = (0:N-1).';

    delta_set(end+1) = make_delta('struct_all_ones', ones(N, 1));
    delta_set(end+1) = make_delta('struct_alternating', (-1) .^ m);
    delta_set(end+1) = make_delta('struct_group_bpsk', group_values(N, M, alternating_group_values(V)));
    delta_set(end+1) = make_delta('struct_halves_same', [ones(N/2, 1); ones(N/2, 1)]);
    delta_set(end+1) = make_delta('struct_halves_opposite', [ones(N/2, 1); -ones(N/2, 1)]);
    delta_set(end+1) = make_delta('struct_group_qpsk', group_values(N, M, cycle_values(V, [1, 1i, -1, -1i])));
    delta_set(end+1) = make_delta('struct_group_qpsk_reverse', group_values(N, M, cycle_values(V, [1, -1i, -1, 1i])));

    [delta_cont, delta_qpsk, delta_bpsk] = build_recursive_delta_from_gps(N, L, d_gps, path_phase);
    delta_set(end+1) = make_delta('recursive_continuous', delta_cont);
    delta_set(end+1) = make_delta('recursive_qpsk', delta_qpsk);
    delta_set(end+1) = make_delta('recursive_bpsk', delta_bpsk);

    oldRngState = rng;
    cleanup = onCleanup(@() rng(oldRngState));
    rng(seed + 1000 * N + 10 * V + L, 'twister');

    if N <= 16
        totalCombos = 2^N;
        combosToUse = min(totalCombos, maxExhaustiveBpsk);
        for combo = 0:combosToUse-1
            bits = dec2bin(combo, N).' - '0';
            values = 2 * bits - 1;
            delta_set(end+1) = make_delta(sprintf('exhaustive_bpsk_%05d', combo + 1), values); %#ok<AGROW>
        end
    else
        for idx = 1:numRandomBpsk
            values = 2 * randi([0, 1], N, 1) - 1;
            delta_set(end+1) = make_delta(sprintf('random_bpsk_%04d', idx), values); %#ok<AGROW>
        end

        qpsk = [1, -1, 1i, -1i].';
        for idx = 1:numRandomQpsk
            values = qpsk(randi(numel(qpsk), N, 1));
            delta_set(end+1) = make_delta(sprintf('random_qpsk_%04d', idx), values); %#ok<AGROW>
        end
    end
end

function item = make_delta(name, delta)
    item.name = name;
    item.delta = delta(:);
    if any(abs(item.delta) == 0)
        error('Delta entries must be nonzero.');
    end
end

function values = group_values(N, M, groupPattern)
    values = zeros(N, 1);
    for group_id = 1:numel(groupPattern)
        values((group_id - 1) * M + 1:group_id * M) = groupPattern(group_id);
    end
end

function values = alternating_group_values(V)
    values = ones(1, V);
    values(2:2:end) = -1;
end

function values = cycle_values(V, baseValues)
    values = zeros(1, V);
    for idx = 1:V
        values(idx) = baseValues(mod(idx - 1, numel(baseValues)) + 1);
    end
end

function value = get_option(options, name, defaultValue)
    if isstruct(options) && isfield(options, name)
        value = options.(name);
    else
        value = defaultValue;
    end
end
