function labels = scheme_display_labels(schemes)
%SCHEME_DISPLAY_LABELS Return plot-friendly scheme names.

    labels = schemes;

    for idx = 1:numel(schemes)
        switch lower(schemes{idx})
            case 'baseline'
                labels{idx} = 'Baseline';
            case 'paper_grouping'
                labels{idx} = 'GPS';
            case 'proposed_grouping'
                labels{idx} = 'Proposed';
            otherwise
                labels{idx} = schemes{idx};
        end
    end
end
