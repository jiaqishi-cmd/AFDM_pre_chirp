function profile = greedy_proposed_profile(symbols, config)
%GREEDY_PROPOSED_PROFILE Select small perturbations by greedy PAPR search.

    N = config.waveform.NumSubcarriers;
    c1 = config.waveform.c1;

    base_profile = proposed_grouping_profile(N, config.pre_chirp);
    selection = greedy_group_papr_selection( ...
        symbols, ...
        N, ...
        c1, ...
        base_profile.candidate_set, ...
        base_profile.group_index);

    profile = base_profile;
    profile.c2 = selection.c2(:);
    profile.selected_candidate_index = selection.selected_candidate_index;
    profile.papr_after_selection = selection.papr;
    profile.selected_signal = selection.signal;
    profile.description = 'Frame-level small-perturbation greedy pre-chirp profile.';
end
