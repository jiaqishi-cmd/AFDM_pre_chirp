function [c2m, d_m] = build_c2m_proposed_pattern(N, V, pattern, c2_base, delta)
%BUILD_C2M_PROPOSED_PATTERN Compatibility wrapper.
%   Prefer afdm.chirp.build_proposed_pattern(N, V, pattern, c2_base, delta)
%   for new code.

    [c2m, d_m] = afdm.chirp.build_proposed_pattern(N, V, pattern, c2_base, delta);
end
