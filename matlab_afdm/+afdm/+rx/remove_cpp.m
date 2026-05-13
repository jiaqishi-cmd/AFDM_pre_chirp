function rData = remove_cpp(rSignal, cppLength)
%REMOVE_CPP Remove the AFDM chirp-periodic prefix.

    if nargin < 2 || isempty(cppLength)
        cppLength = 0;
    end
    if cppLength < 0 || cppLength > numel(rSignal)
        error('CPP length must be between 0 and the signal length.');
    end

    rData = rSignal(cppLength + 1 : end);
end
