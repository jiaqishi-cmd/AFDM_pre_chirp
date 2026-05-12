function delta_set = generate_delta_set_for_dmin(N, M, options)
%GENERATE_DELTA_SET_FOR_DMIN 生成用于有限星座最小距离搜索的 delta 候选。
%   delta 只作为差分方向搜索，主实验使用结构化方向加随机 BPSK/QPSK 相位方向。

    if nargin < 3
        options = struct();
    end
    numBpsk = get_option(options, 'num_random_bpsk', 2000);
    numQpsk = get_option(options, 'num_random_qpsk', 2000);
    seed = get_option(options, 'seed', 20260508);

    oldRngState = rng;
    cleanup = onCleanup(@() rng(oldRngState)); %#ok<NASGU>
    rng(seed, 'twister');

    delta_mat = [];
    delta_type = strings(0, 1);
    delta_index = zeros(0, 1);

    append_delta('all_ones', ones(N, 1), 1);
    append_delta('group_alt', group_delta(N, M, [1, -1, 1, -1]), 1);
    append_delta('half_signs', group_delta(N, M, [1, 1, -1, -1]), 1);
    append_delta('alternating', (-1) .^ (0:N-1).', 1);
    append_delta('qpsk_group_1', group_delta(N, M, [1, 1i, -1, -1i]), 1);
    append_delta('qpsk_group_2', group_delta(N, M, [1, -1i, -1, 1i]), 1);

    for idx = 1:numBpsk
        append_delta('random_bpsk', 2 * randi([0, 1], N, 1) - 1, idx);
    end

    alphabet = [1, -1, 1i, -1i].';
    for idx = 1:numQpsk
        append_delta('random_qpsk', alphabet(randi(numel(alphabet), N, 1)), idx);
    end

    delta_set.delta = delta_mat;
    delta_set.type = delta_type;
    delta_set.index = delta_index;
    delta_set.norm2 = sum(abs(delta_mat).^2, 1);

    function append_delta(name, delta, idx)
        delta_mat(:, end+1) = delta(:); %#ok<AGROW>
        delta_type(end+1, 1) = string(name); %#ok<AGROW>
        delta_index(end+1, 1) = idx; %#ok<AGROW>
    end
end

function delta = group_delta(N, M, values)
    delta = zeros(N, 1);
    numGroups = N / M;
    for groupIdx = 1:numGroups
        value = values(mod(groupIdx - 1, numel(values)) + 1);
        firstIdx = (groupIdx - 1) * M + 1;
        lastIdx = groupIdx * M;
        delta(firstIdx:lastIdx) = value;
    end
end

function value = get_option(options, name, defaultValue)
    if isstruct(options) && isfield(options, name)
        value = options.(name);
    else
        value = defaultValue;
    end
end
