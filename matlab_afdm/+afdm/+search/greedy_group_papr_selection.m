function selection = greedy_group_papr_selection(symbols, numSubcarriers, c1, candidateSet, groupIndex, oversamplingFactor)
%GREEDY_GROUP_PAPR_SELECTION Greedily choose group candidates to reduce PAPR.

    if nargin < 6 || isempty(oversamplingFactor)
        oversamplingFactor = 1;
    end

    numGroups = max(groupIndex);
    numCandidates = size(candidateSet, 2);

    c2Vec = candidateSet(:, 1);
    bestSignal = afdm.search.full_waveform(symbols, c2Vec, oversamplingFactor);
    bestPapr = afdm.tx.compute_papr(bestSignal);
    selectedCandidateIndex = ones(numGroups, 1);

    for groupId = 1:numGroups
        indices = find(groupIndex == groupId);
        groupBestC2 = c2Vec(indices);
        groupBestPapr = bestPapr;
        groupBestSignal = bestSignal;
        groupBestCandidate = selectedCandidateIndex(groupId);

        for candidateId = 2:numCandidates
            c2Trial = c2Vec;
            c2Trial(indices) = candidateSet(indices, candidateId);

            signalTrial = afdm.search.full_waveform(symbols, c2Trial, oversamplingFactor);
            paprTrial = afdm.tx.compute_papr(signalTrial);

            if paprTrial < groupBestPapr
                groupBestC2 = c2Trial(indices);
                groupBestPapr = paprTrial;
                groupBestSignal = signalTrial;
                groupBestCandidate = candidateId;
            end
        end

        c2Vec(indices) = groupBestC2;
        bestPapr = groupBestPapr;
        bestSignal = groupBestSignal;
        selectedCandidateIndex(groupId) = groupBestCandidate;
    end

    selection.c2 = c2Vec;
    selection.papr = bestPapr;
    selection.signal = bestSignal;
    selection.selected_candidate_index = selectedCandidateIndex;
    selection.oversampling_factor = oversamplingFactor;
end
