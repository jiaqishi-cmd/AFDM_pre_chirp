function config = apply_pre_chirp_scheme(config, scheme)
%APPLY_PRE_CHIRP_SCHEME Update config with a selected pre-chirp profile.

    if nargin < 2 || isempty(scheme)
        scheme = config.pre_chirp.scheme;
    end

    config.pre_chirp.scheme = scheme;
    config.pre_chirp.profile = generate_pre_chirp_profile( ...
        scheme, ...
        config.waveform.NumSubcarriers, ...
        config.pre_chirp);
    config.waveform.c2 = config.pre_chirp.profile.c2;
end
