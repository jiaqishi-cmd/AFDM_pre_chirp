function [signal_cpp, papr, bits, tx_state] = afdm_tx_engine(config)
%AFDM_TX_ENGINE Compatibility wrapper.
%   Prefer afdm.tx.engine(config) for new code.

    if nargin < 1
        config = [];
    end
    [signal_cpp, papr, bits, tx_state] = afdm.tx.engine(config);
end
