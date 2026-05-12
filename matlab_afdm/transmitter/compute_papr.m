function papr = compute_papr(signal)
%COMPUTE_PAPR Compatibility wrapper.
%   Prefer afdm.tx.compute_papr(signal) for new code.

    papr = afdm.tx.compute_papr(signal);
end
