function signal_cpp = add_cpp(signal, cppLength, c1)
%ADD_CPP Compatibility wrapper.
%   Prefer afdm.tx.add_cpp(signal, cppLength, c1).

    if nargin < 3
        c1 = [];
    end
    signal_cpp = afdm.tx.add_cpp(signal, cppLength, c1);
end
