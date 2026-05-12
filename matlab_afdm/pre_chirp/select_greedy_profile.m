function profile = select_greedy_profile(symbols, config, scheme)
%SELECT_GREEDY_PROFILE Apply grouped greedy PAPR selection for a scheme.

    N = config.waveform.NumSubcarriers;
    c1 = config.waveform.c1;

    base_profile = afdm.chirp.build_profile(scheme, N, config.pre_chirp);
    selection = greedy_group_papr_selection( ...
        symbols, ...
        N, ...
        c1, ...
        base_profile.candidate_set, ...
        base_profile.group_index);

    profile = base_profile;
    profile.c2 = selection.c2(:);
    profile.selection.selected_candidate_index = selection.selected_candidate_index;
    profile.selection.papr = selection.papr;
    profile.selection.signal = selection.signal;
    profile.selected_candidate_index = selection.selected_candidate_index;
    profile.papr_after_selection = selection.papr;
    profile.selected_signal = selection.signal;
    profile.description = sprintf('Frame-level %s greedy pre-chirp profile.', scheme);
end
