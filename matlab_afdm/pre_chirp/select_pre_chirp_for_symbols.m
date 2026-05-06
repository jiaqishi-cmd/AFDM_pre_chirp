function profile = select_pre_chirp_for_symbols(symbols, config)
%SELECT_PRE_CHIRP_FOR_SYMBOLS Select the frame-level pre-chirp profile.

    scheme = config.pre_chirp.scheme;

    switch lower(scheme)
        case 'baseline'
            profile = config.pre_chirp.profile;
        case 'paper_grouping'
            profile = select_greedy_profile(symbols, config, scheme);
        case 'proposed_grouping'
            profile = select_greedy_profile(symbols, config, scheme);
        otherwise
            error('Unsupported pre-chirp scheme: %s', scheme);
    end
end
