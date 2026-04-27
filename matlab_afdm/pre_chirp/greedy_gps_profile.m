function profile = greedy_gps_profile(symbols, config)
%GREEDY_GPS_PROFILE Select c2 by the GPS greedy PAPR minimization rule.

    N = config.waveform.NumSubcarriers;
    c1 = config.waveform.c1;

    base_profile = paper_grouping_profile(N, config.pre_chirp);
    candidate_set = base_profile.candidate_set;
    group_index = base_profile.group_index;
    num_groups = max(group_index);
    num_candidates = size(candidate_set, 2);

    c2_vec = candidate_set(:, 1);
    best_signal = idaft_mod(symbols, N, c1, c2_vec);
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

            signal_trial = idaft_mod(symbols, N, c1, c2_trial);
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

    profile = base_profile;
    profile.c2 = c2_vec(:);
    profile.selected_candidate_index = selected_candidate_index;
    profile.papr_after_selection = best_papr;
    profile.selected_signal = best_signal;
    profile.description = 'Frame-level GPS greedy pre-chirp profile.';
end
