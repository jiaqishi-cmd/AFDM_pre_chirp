function [err_bits, total_bits] = compute_bit_errors(x_dec, x_ref)
%COMPUTE_BIT_ERRORS Count bit mismatches between decisions and references.
%   [err_bits, total_bits] = compute_bit_errors(x_dec, x_ref)

    if numel(x_dec) ~= numel(x_ref)
        error('Decision bits and reference bits must have the same length.');
    end

    err_bits = sum(x_dec(:) ~= x_ref(:));
    total_bits = numel(x_ref);
end
