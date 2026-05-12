function [c2m, d_m] = build_c2m_gps_pattern(N, V, pattern)
%BUILD_C2M_GPS_PATTERN Compatibility wrapper.
%   Prefer afdm.chirp.build_gps_pattern(N, V, pattern) for new code.

    [c2m, d_m] = afdm.chirp.build_gps_pattern(N, V, pattern);
end
