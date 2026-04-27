function profile = greedy_gps_profile(symbols, config)
%GREEDY_GPS_PROFILE Select c2 by the GPS greedy PAPR minimization rule.

    N = config.waveform.NumSubcarriers;
    c1 = config.waveform.c1;

    base_profile = paper_grouping_profile(N, config.pre_chirp);
    candidate_set = base_profile.candidate_set;
    group_index = base_profile.group_index;
    num_groups = max(group_index);
    num_candidates = size(candidate_set, 2);

    selection = greedy_group_papr_selection(symbols, N, c1, candidate_set, group_index);

    profile = base_profile;
    profile.c2 = selection.c2(:);
    profile.selected_candidate_index = selection.selected_candidate_index;
    profile.papr_after_selection = selection.papr;
    profile.selected_signal = selection.signal;
    profile.description = 'Frame-level GPS greedy pre-chirp profile.';
end
