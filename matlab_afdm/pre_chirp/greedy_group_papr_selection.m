function selection = greedy_group_papr_selection(symbols, numSubcarriers, c1, candidate_set, group_index)
%GREEDY_GROUP_PAPR_SELECTION Greedily choose group candidates to reduce PAPR.

    num_groups = max(group_index);
    num_candidates = size(candidate_set, 2);

    c2_vec = candidate_set(:, 1);
    best_signal = idaft_mod(symbols, numSubcarriers, c1, c2_vec);
    best_papr = compute_papr(best_signal);
    selected_candidate_index = ones(num_groups, 1);

    for group_id = 1:num_groups
        indices = find(group_index == group_id);
        group_best_c2 = c2_vec(indices);
        group_best_papr = best_papr;
        group_best_signal = best_signal;
        group_best_candidate = selected_candidate_index(group_id);

        for candidate_id = 2:num_candidates
            c2_trial = c2_vec;
            c2_trial(indices) = candidate_set(indices, candidate_id);

            signal_trial = idaft_mod(symbols, numSubcarriers, c1, c2_trial);
            papr_trial = compute_papr(signal_trial);

            if papr_trial < group_best_papr
                group_best_c2 = c2_trial(indices);
                group_best_papr = papr_trial;
                group_best_signal = signal_trial;
                group_best_candidate = candidate_id;
            end
        end

        c2_vec(indices) = group_best_c2;
        best_papr = group_best_papr;
        best_signal = group_best_signal;
        selected_candidate_index(group_id) = group_best_candidate;
    end

    selection.c2 = c2_vec;
    selection.papr = best_papr;
    selection.signal = best_signal;
    selection.selected_candidate_index = selected_candidate_index;
end
