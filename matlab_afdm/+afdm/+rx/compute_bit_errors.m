function [errBits, totalBits] = compute_bit_errors(xDec, xRef)
%COMPUTE_BIT_ERRORS Count bit mismatches between decisions and references.

    if numel(xDec) ~= numel(xRef)
        error('Decision bits and reference bits must have the same length.');
    end

    errBits = sum(xDec(:) ~= xRef(:));
    totalBits = numel(xRef);
end
