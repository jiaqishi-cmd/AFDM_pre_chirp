function out = reuse_beam_search(symbols, baseC2, offsets, groupIndex, searchOs, finalOs, beamWidth, topK)
%REUSE_BEAM_SEARCH Beam search using partial waveform reuse.

    V = max(groupIndex);
    W = numel(offsets);
    M = numel(symbols);

    sPartSearch = afdm.search.precompute_partial_waveforms(symbols, baseC2, offsets, groupIndex, searchOs);
    numIfft = V * W;
    initPattern = ceil(W / 2) * ones(1, V);
    initWaveform = afdm.search.combine_partial_waveform(sPartSearch, initPattern);
    states = make_initial_state(M, V, ceil(W / 2));
    states.waveform = initWaveform;
    states.metric = afdm.tx.compute_papr(initWaveform);
    evalSearch = 0;

    for groupId = 1:V
        expanded = repmat(make_state(M, V), 1, numel(states) * W);
        outIdx = 0;
        for stateIdx = 1:numel(states)
            oldCand = states(stateIdx).pattern(groupId);
            for candId = 1:W
                outIdx = outIdx + 1;
                pattern = states(stateIdx).pattern;
                pattern(groupId) = candId;
                waveform = states(stateIdx).waveform - sPartSearch{groupId, oldCand} + sPartSearch{groupId, candId};
                expanded(outIdx).pattern = pattern;
                expanded(outIdx).waveform = waveform;
                expanded(outIdx).metric = afdm.tx.compute_papr(waveform);
                evalSearch = evalSearch + 1;
            end
        end
        [~, order] = sort([expanded.metric], 'ascend');
        states = expanded(order(1:min(beamWidth, numel(order))));
    end

    sPartFinal = afdm.search.precompute_partial_waveforms(symbols, baseC2, offsets, groupIndex, finalOs);
    numIfft = numIfft + V * W;
    finalK = min(topK, numel(states));
    finalPapr = zeros(1, finalK);
    for idx = 1:finalK
        s = afdm.search.combine_partial_waveform(sPartFinal, states(idx).pattern);
        finalPapr(idx) = afdm.tx.compute_papr(s);
    end

    [bestPapr, bestIdx] = min(finalPapr);
    out.papr = bestPapr;
    out.pattern = states(bestIdx).pattern;
    out.eval_search = evalSearch;
    out.eval_final = finalK;
    out.eval_count = evalSearch + finalK;
    out.num_ifft = numIfft;
end

function state = make_initial_state(M, V, candId)
    state = make_state(M, V);
    state.pattern = candId * ones(1, V);
    state.metric = Inf;
end

function state = make_state(M, V)
    state = struct('pattern', zeros(1, V), 'metric', Inf, 'waveform', zeros(M, 1));
end
