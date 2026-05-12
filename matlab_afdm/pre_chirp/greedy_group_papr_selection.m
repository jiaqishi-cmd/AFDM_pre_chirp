function selection = greedy_group_papr_selection(symbols, numSubcarriers, c1, candidate_set, group_index)
%GREEDY_GROUP_PAPR_SELECTION Compatibility wrapper.
%   Prefer afdm.search.greedy_group_papr_selection for new code.

    selection = afdm.search.greedy_group_papr_selection( ...
        symbols, numSubcarriers, c1, candidate_set, group_index);
end
